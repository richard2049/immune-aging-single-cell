from __future__ import annotations

import argparse
import json
import re
import zlib
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Any

import numpy as np
import pandas as pd

from .utils import load_config

REQUIRED_RESULT_COLUMNS = {
    "gene_id",
    "cell_type",
    "analysis_tier",
    "model_scope",
    "model_formula",
    "log2_fc_per_10_years",
    "average_log_cpm",
    "ql_f_statistic",
    "p_value",
    "fdr_within_celltype",
    "fdr_global",
    "n_replicates",
    "age_span_years",
}
LONGITUDINAL_RESULT_COLUMNS = {
    "gene_id",
    "cell_type",
    "analysis_tier",
    "model_scope",
    "model_formula",
    "log2_fc_per_10_years",
    "average_log_cpm",
    "moderated_t_statistic",
    "p_value",
    "fdr_within_celltype",
    "fdr_global",
    "n_sample_units",
    "n_subjects",
    "age_span_years",
}
REQUIRED_PROFILE_COLUMNS = {
    "profile_id",
    "biological_replicate_id",
    "cell_type",
    "age",
    "age_decade",
    "sex",
    "batch",
    "analysis_tier",
    "n_cells",
    "library_size",
}


def _bh_adjust(values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(values, dtype=float)
    if p_values.ndim != 1 or not np.isfinite(p_values).all():
        raise ValueError("BH adjustment requires a finite one-dimensional array.")
    if ((p_values < 0) | (p_values > 1)).any():
        raise ValueError("P-values must be between zero and one.")

    order = np.argsort(p_values, kind="stable")
    ranked = p_values[order]
    ranks = np.arange(1, len(ranked) + 1, dtype=float)
    adjusted = ranked * len(ranked) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _resolve_manifest_path(output_dir: Path, value: str) -> Path:
    raw = str(value).strip()
    windows_path = PureWindowsPath(raw)
    if not raw or raw.startswith(("/", "\\")) or windows_path.drive or Path(raw).is_absolute():
        raise ValueError(f"Manifest path must be relative: {value!r}")

    candidate = (output_dir / Path(raw)).resolve()
    root = output_dir.resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ValueError(
            f"Manifest path escapes the pseudobulk output directory: {value!r}"
        ) from error
    return candidate


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _inspect_r_graphics_pdf(path: Path, expected_pages: int) -> tuple[bool, str]:
    content = path.read_bytes()
    if not content.startswith(b"%PDF-") or b"%%EOF" not in content[-1024:]:
        return False, "The diagnostic file lacks a complete PDF header or trailer."

    page_count = len(re.findall(rb"/Type\s*/Page\b", content))
    errors: list[str] = []
    stream_count = 0
    stream_pattern = re.compile(rb"<<(.*?)>>\s*stream\r?\n", re.DOTALL)
    for match in stream_pattern.finditer(content):
        dictionary = match.group(1)
        if b"/FlateDecode" not in dictionary:
            continue
        length_match = re.search(rb"/Length\s+(\d+)\b", dictionary)
        if length_match is None:
            errors.append("compressed stream with non-literal length")
            continue
        stream_count += 1
        stream_length = int(length_match.group(1))
        stream = content[match.end() : match.end() + stream_length]
        if len(stream) != stream_length:
            errors.append(f"truncated compressed stream {stream_count}")
            continue
        try:
            zlib.decompress(stream)
        except zlib.error as error:
            errors.append(f"invalid compressed stream {stream_count}: {error}")

    passed = page_count == expected_pages and stream_count >= expected_pages and not errors
    detail = (
        f"Validated {page_count} PDF pages and {stream_count} compressed streams."
        if passed
        else (
            f"Expected {expected_pages} pages; observed {page_count} pages and "
            f"{stream_count} compressed streams. " + "; ".join(errors)
        ).strip()
    )
    return passed, detail


def _eligible_analysis_sets(
    contract: dict[str, Any],
    eligibility: pd.DataFrame,
) -> tuple[set[str], set[str], bool]:
    configured_tiers = {
        **{str(cell_type): "primary" for cell_type in contract["primary_cell_types"]},
        **{str(cell_type): "exploratory" for cell_type in contract["exploratory_cell_types"]},
    }
    required = {"cell_type", "analysis_tier", "aggregation_status", "exclusion_reasons"}
    if not required.issubset(eligibility.columns):
        return set(), set(), False

    cell_types = eligibility["cell_type"].astype(str)
    statuses = eligibility["aggregation_status"].astype(str)
    reasons = eligibility["exclusion_reasons"].fillna("").astype(str).str.strip()
    observed_tiers = dict(zip(cell_types, eligibility["analysis_tier"].astype(str)))
    contract_ok = (
        cell_types.is_unique
        and set(cell_types) == set(configured_tiers)
        and statuses.isin({"included", "excluded"}).all()
        and observed_tiers == configured_tiers
        and reasons.loc[statuses.eq("excluded")].ne("").all()
    )
    included = set(cell_types.loc[statuses.eq("included")])
    eligible_primary = {
        cell_type for cell_type in included if configured_tiers.get(cell_type) == "primary"
    }
    eligible_exploratory = {
        cell_type for cell_type in included if configured_tiers.get(cell_type) == "exploratory"
    }
    return eligible_primary, eligible_exploratory, bool(contract_ok)


def validate(
    config_path: str | Path,
    output_dir: str | Path,
    report_path: str | Path,
) -> dict[str, Any]:
    config_path = Path(config_path)
    output_dir = Path(output_dir)
    report_path = Path(report_path)
    config = load_config(config_path)
    contract = config["pseudobulk_de"]
    runtime_contract = contract["runtime"]
    longitudinal = bool(contract.get("subject_col"))

    paths = {
        "aggregation_audit": output_dir / "aggregation_audit.json",
        "counts": output_dir / "pseudobulk_counts.mtx.gz",
        "profiles": output_dir / "pseudobulk_profiles.csv",
        "genes": output_dir / "pseudobulk_genes.csv",
        "eligibility": output_dir / "celltype_eligibility.csv",
        "results": output_dir / "gene_level_results.csv.gz",
        "manifest": output_dir / "celltype_result_manifest.csv",
        "diagnostics": output_dir / "design_diagnostics.csv",
        "diagnostic_plots": output_dir
        / ("dream_diagnostic_plots.pdf" if longitudinal else "edgeR_diagnostic_plots.pdf"),
        "runtime_versions": output_dir / "runtime_versions.csv",
        "session_info": output_dir / "session_info.txt",
        "model_log": output_dir / ("dream.log" if longitudinal else "edgeR.log"),
    }
    if longitudinal:
        paths["sensitivity_results"] = output_dir / "one_sample_per_subject_results.csv.gz"
    checks: list[dict[str, Any]] = []

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    missing = [
        str(path) for path in paths.values() if not path.is_file() or path.stat().st_size == 0
    ]
    record(
        "required_outputs_exist",
        not missing,
        "All required outputs are present and non-empty."
        if not missing
        else "Missing or empty outputs: " + ", ".join(missing),
    )
    if missing:
        report = {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "config_path": str(config_path),
            "output_dir": str(output_dir),
            "passed": False,
            "requires_manual_review": True,
            "requires_technical_manual_review": True,
            "biological_interpretation_required": True,
            "interpretation_status": ("technical_validation_only_no_biological_claims"),
            "checks": checks,
            "summary": {},
        }
        _write_json(report_path, report)
        raise ValueError("Pseudobulk technical validation failed: required outputs are missing.")

    audit = json.loads(paths["aggregation_audit"].read_text(encoding="utf-8"))
    profiles = pd.read_csv(paths["profiles"])
    genes = pd.read_csv(paths["genes"])
    eligibility = pd.read_csv(paths["eligibility"])
    results = pd.read_csv(paths["results"])
    sensitivity_results = (
        pd.read_csv(paths["sensitivity_results"]) if longitudinal else pd.DataFrame()
    )
    manifest = pd.read_csv(paths["manifest"])
    diagnostics = pd.read_csv(paths["diagnostics"])
    versions = pd.read_csv(paths["runtime_versions"], dtype=str).fillna("")

    log_errors: list[str] = []
    log_contracts = [(paths["model_log"], "return_code: 0")]
    aggregation_log = output_dir / "aggregation.log"
    if aggregation_log.exists():
        log_contracts.append((aggregation_log, "exit_code: 0"))
    for log_path, completion_marker in log_contracts:
        try:
            log_text = log_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            log_errors.append(f"Log is not valid UTF-8: {log_path}")
            continue
        if completion_marker not in log_text:
            log_errors.append(f"Completion marker {completion_marker!r} is absent: {log_path}")
    record(
        "execution_logs",
        not log_errors,
        "Execution logs are UTF-8 and contain successful completion markers."
        if not log_errors
        else "; ".join(log_errors),
    )

    record(
        "aggregation_contract",
        bool(audit.get("count_conservation_passed"))
        and not bool(audit.get("dense_cell_by_gene_conversion"))
        and audit.get("matrix_orientation") == "profiles_by_genes"
        and audit.get("n_profiles") == len(profiles)
        and audit.get("n_genes") == len(genes),
        (
            "Aggregation reports exact count conservation, sparse processing, "
            "and matching dimensions."
        ),
    )

    required_profile_columns = (
        (REQUIRED_PROFILE_COLUMNS - {"biological_replicate_id"}) | {"sample_unit_id", "subject_id"}
        if longitudinal
        else REQUIRED_PROFILE_COLUMNS
    )
    profile_id_col = "sample_unit_id" if longitudinal else "biological_replicate_id"
    profile_columns_ok = required_profile_columns.issubset(profiles.columns)
    profile_keys_ok = (
        profile_columns_ok
        and profiles["profile_id"].notna().all()
        and not profiles["profile_id"].duplicated().any()
        and not profiles[[profile_id_col, "cell_type"]].duplicated().any()
        and profiles[list(required_profile_columns)].notna().all().all()
        and profiles["n_cells"].ge(contract["min_cells_per_pseudobulk"]).all()
    )
    record(
        "profile_metadata_integrity",
        profile_keys_ok,
        f"Checked {len(profiles):,} unique replicate-by-cell-type profiles.",
    )

    configured_primary = set(contract["primary_cell_types"])
    configured_exploratory = set(contract["exploratory_cell_types"])
    primary, exploratory, eligibility_contract_ok = _eligible_analysis_sets(contract, eligibility)
    approved = primary | exploratory
    observed_tiers = {
        tier: set(group["cell_type"])
        for tier, group in profiles.groupby("analysis_tier", sort=False)
    }
    tier_contract_ok = (
        observed_tiers.get("primary", set()) == primary
        and observed_tiers.get("exploratory", set()) == exploratory
        and set(profiles["cell_type"]) == approved
    )
    record(
        "analysis_tiers_match_contract",
        tier_contract_ok,
        (
            f"Observed {len(primary)} eligible primary and "
            f"{len(exploratory)} eligible exploratory cell types."
        ),
    )

    included = eligibility.loc[eligibility["aggregation_status"].eq("included"), "cell_type"]
    record(
        "eligibility_matches_contract",
        eligibility_contract_ok and set(included) == approved and not included.duplicated().any(),
        (
            f"Eligibility records {len(configured_primary)} configured primary and "
            f"{len(configured_exploratory)} configured exploratory cell types; "
            f"{len(included)} are approved for pseudobulk modeling."
        ),
    )

    required_result_columns = (
        LONGITUDINAL_RESULT_COLUMNS if longitudinal else REQUIRED_RESULT_COLUMNS
    )
    result_columns_ok = required_result_columns.issubset(results.columns)
    numeric_columns = [
        "log2_fc_per_10_years",
        "average_log_cpm",
        "moderated_t_statistic" if longitudinal else "ql_f_statistic",
        "p_value",
        "fdr_within_celltype",
        "fdr_global",
        "n_sample_units" if longitudinal else "n_replicates",
        "n_subjects" if longitudinal else "n_replicates",
        "age_span_years",
    ]
    numeric_ok = (
        result_columns_ok and np.isfinite(results[numeric_columns].to_numpy(dtype=float)).all()
    )
    probabilities_ok = result_columns_ok and (
        results[["p_value", "fdr_within_celltype", "fdr_global"]].ge(0).all().all()
        and results[["p_value", "fdr_within_celltype", "fdr_global"]].le(1).all().all()
    )
    unique_tests = result_columns_ok and not results[["gene_id", "cell_type"]].duplicated().any()
    record(
        "gene_level_result_integrity",
        result_columns_ok and numeric_ok and probabilities_ok and unique_tests,
        f"Checked {len(results):,} unique gene-by-cell-type tests.",
    )

    sensitivity_ok = True
    if longitudinal:
        sensitivity_required = LONGITUDINAL_RESULT_COLUMNS | {"selection_rule"}
        sensitivity_columns_ok = sensitivity_required.issubset(sensitivity_results.columns)
        sensitivity_numeric_ok = (
            sensitivity_columns_ok
            and np.isfinite(sensitivity_results[numeric_columns].to_numpy(dtype=float)).all()
        )
        sensitivity_probabilities_ok = sensitivity_columns_ok and (
            sensitivity_results[["p_value", "fdr_within_celltype", "fdr_global"]].ge(0).all().all()
            and sensitivity_results[["p_value", "fdr_within_celltype", "fdr_global"]]
            .le(1)
            .all()
            .all()
        )
        sensitivity_unique = (
            sensitivity_columns_ok
            and not sensitivity_results[["gene_id", "cell_type"]].duplicated().any()
        )
        sensitivity_ok = (
            sensitivity_columns_ok
            and sensitivity_results[list(sensitivity_required)].notna().all().all()
            and sensitivity_numeric_ok
            and sensitivity_probabilities_ok
            and sensitivity_unique
            and set(sensitivity_results["cell_type"]) == approved
            and sensitivity_results["model_scope"].eq("one_sample_per_subject").all()
            and sensitivity_results["selection_rule"].eq(contract["sensitivity_sample_rule"]).all()
            and sensitivity_results["n_sample_units"].eq(sensitivity_results["n_subjects"]).all()
            and sensitivity_results["model_formula"].str.contains("age_decade", regex=False).all()
            and not sensitivity_results["model_formula"].str.contains("|", regex=False).any()
        )
        record(
            "one_sample_per_subject_sensitivity",
            sensitivity_ok,
            (
                f"Checked {len(sensitivity_results):,} deterministic sensitivity tests "
                f"using {contract['sensitivity_sample_rule']}."
            ),
        )

    expected_tiers = {
        **{cell_type: "primary" for cell_type in primary},
        **{cell_type: "exploratory" for cell_type in exploratory},
    }
    result_tiers_ok = all(
        set(group["analysis_tier"]) == {expected_tiers[cell_type]}
        for cell_type, group in results.groupby("cell_type", sort=False)
    )
    accepted_scope = "adjusted_repeated_measures" if longitudinal else "adjusted_primary"
    model_formula_ok = (
        results.loc[results["model_scope"].eq(accepted_scope), "model_formula"]
        .eq(contract["primary_formula"])
        .all()
    )
    record(
        "result_analysis_contract",
        set(results["cell_type"]) == approved and result_tiers_ok and model_formula_ok,
        "Result tiers and adjusted formulas match the approved contract.",
    )

    diagnostic_cell_types = set(diagnostics["cell_type"])
    expected_diagnostic_scopes = (
        {accepted_scope, "one_sample_per_subject"} if longitudinal else {accepted_scope}
    )
    primary_diagnostics = diagnostics[diagnostics["model_scope"].eq(accepted_scope)]
    sensitivity_diagnostics = diagnostics[diagnostics["model_scope"].eq("one_sample_per_subject")]
    diagnostics_ok = (
        diagnostic_cell_types == approved
        and set(diagnostics["model_scope"]) == expected_diagnostic_scopes
        and not diagnostics[["cell_type", "model_scope"]].duplicated().any()
        and diagnostics["status"].eq("completed").all()
        and diagnostics["design_rank"].eq(diagnostics["design_columns"]).all()
        and diagnostics["residual_df"].ge(contract["min_residual_df"]).all()
        and primary_diagnostics["n_genes_tested"].sum() == len(results)
        and (
            not longitudinal
            or sensitivity_diagnostics["n_genes_tested"].sum() == len(sensitivity_results)
        )
    )
    record(
        "model_diagnostics",
        diagnostics_ok,
        (
            "All approved cell types completed with full-rank designs and "
            "sufficient residual degrees of freedom."
        ),
    )

    pdf_ok, pdf_detail = _inspect_r_graphics_pdf(
        paths["diagnostic_plots"],
        expected_pages=len(primary_diagnostics),
    )
    record("diagnostic_plot_pdf_integrity", pdf_ok, pdf_detail)

    fallback = diagnostics.loc[
        ~diagnostics["model_scope"].isin(expected_diagnostic_scopes), "cell_type"
    ].tolist()
    record(
        "approved_primary_models",
        not fallback,
        "All cell types used the approved adjusted model."
        if not fallback
        else "Manual review required for non-primary model scopes: " + ", ".join(fallback),
    )

    global_recomputed = _bh_adjust(results["p_value"].to_numpy())
    global_difference = float(np.max(np.abs(global_recomputed - results["fdr_global"].to_numpy())))
    within_difference = 0.0
    for _, indices in results.groupby("cell_type", sort=False).groups.items():
        index = np.asarray(list(indices), dtype=int)
        expected = _bh_adjust(results.loc[index, "p_value"].to_numpy())
        observed = results.loc[index, "fdr_within_celltype"].to_numpy()
        within_difference = max(
            within_difference,
            float(np.max(np.abs(expected - observed))),
        )
    fdr_tolerance = 1e-12
    record(
        "multiple_testing_recalculation",
        global_difference <= fdr_tolerance and within_difference <= fdr_tolerance,
        (
            f"Maximum absolute differences: global={global_difference:.3g}, "
            f"within-cell-type={within_difference:.3g}."
        ),
    )

    manifest_ok = (
        set(manifest["cell_type"]) == approved
        and manifest["cell_type"].is_unique
        and int(manifest["n_genes_tested"].sum()) == len(results)
    )
    manifest_errors: list[str] = []
    if manifest_ok:
        for row in manifest.itertuples(index=False):
            try:
                result_path = _resolve_manifest_path(output_dir, row.result_path)
            except ValueError as error:
                manifest_errors.append(str(error))
                continue
            if not result_path.is_file() or result_path.stat().st_size == 0:
                manifest_errors.append(f"Missing result table: {result_path}")
                continue
            table = pd.read_csv(result_path)
            if len(table) != int(row.n_genes_tested):
                manifest_errors.append(f"Row-count mismatch: {result_path}")
            if set(table["cell_type"]) != {row.cell_type}:
                manifest_errors.append(f"Cell-type mismatch: {result_path}")
    record(
        "cell_type_manifest",
        manifest_ok and not manifest_errors,
        (
            "All manifest paths are portable and all indexed tables match "
            "their declared cell type and row count."
        )
        if manifest_ok and not manifest_errors
        else "; ".join(manifest_errors) or "Manifest rows do not match the result contract.",
    )

    observed_versions = dict(zip(versions["component"], versions["observed"]))
    versions_ok = (
        observed_versions.get("R", "").startswith(runtime_contract["expected_r_version"])
        and observed_versions.get("Bioconductor")
        == runtime_contract["expected_bioconductor_version"]
        and observed_versions.get("edgeR") == runtime_contract["expected_edger_version"]
        and (
            not longitudinal
            or observed_versions.get("variancePartition")
            == runtime_contract["expected_variance_partition_version"]
        )
    )
    record(
        "runtime_versions",
        versions_ok,
        (
            f"Observed R {observed_versions.get('R')}, Bioconductor "
            f"{observed_versions.get('Bioconductor')}, edgeR "
            f"{observed_versions.get('edgeR')}, and variancePartition "
            f"{observed_versions.get('variancePartition')}."
        ),
    )

    temporary_outputs = [str(path) for path in output_dir.rglob("*.tmp") if path != report_path]
    record(
        "no_incomplete_outputs",
        not temporary_outputs,
        "No temporary output files remain."
        if not temporary_outputs
        else "Temporary outputs remain: " + ", ".join(temporary_outputs),
    )

    checks_passed = all(check["passed"] for check in checks)
    summary = {
        "n_profiles": int(len(profiles)),
        "n_sample_units": int(profiles[profile_id_col].nunique()),
        "n_subjects": int(
            profiles["subject_id"].nunique() if longitudinal else profiles[profile_id_col].nunique()
        ),
        "n_genes_input": int(len(genes)),
        "n_cell_types_tested": int(results["cell_type"].nunique()),
        "n_configured_primary_cell_types": int(len(configured_primary)),
        "n_configured_exploratory_cell_types": int(len(configured_exploratory)),
        "n_primary_cell_types": int(len(primary)),
        "n_exploratory_cell_types": int(len(exploratory)),
        "n_gene_celltype_tests": int(len(results)),
        "n_sensitivity_gene_celltype_tests": int(len(sensitivity_results)),
        "n_global_fdr_lt_0_05": int(results["fdr_global"].lt(0.05).sum()),
        "n_within_celltype_fdr_lt_0_05": int(results["fdr_within_celltype"].lt(0.05).sum()),
        "model_scopes": {
            str(key): int(value) for key, value in diagnostics["model_scope"].value_counts().items()
        },
        "tests_by_cell_type": {
            str(key): int(value) for key, value in results.groupby("cell_type").size().items()
        },
    }
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path),
        "output_dir": str(output_dir),
        "passed": checks_passed,
        "requires_manual_review": True,
        "requires_technical_manual_review": bool(fallback) or not checks_passed,
        "biological_interpretation_required": True,
        "interpretation_status": ("technical_validation_only_no_biological_claims"),
        "checks": checks,
        "summary": summary,
    }
    _write_json(report_path, report)
    if not checks_passed:
        failed = [check["name"] for check in checks if not check["passed"]]
        raise ValueError("Pseudobulk technical validation failed: " + ", ".join(failed))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description=("Validate donor-aware pseudobulk DE outputs without interpreting biology.")
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    report = validate(args.config, args.outdir, args.report)
    summary = report["summary"]
    print(
        "[validate_pseudobulk_de] passed: "
        f"{summary['n_cell_types_tested']} cell types, "
        f"{summary['n_gene_celltype_tests']:,} tests"
    )


if __name__ == "__main__":
    main()
