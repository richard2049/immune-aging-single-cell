from __future__ import annotations

import argparse
import json
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

    paths = {
        "aggregation_audit": output_dir / "aggregation_audit.json",
        "counts": output_dir / "pseudobulk_counts.mtx.gz",
        "profiles": output_dir / "pseudobulk_profiles.csv",
        "genes": output_dir / "pseudobulk_genes.csv",
        "eligibility": output_dir / "celltype_eligibility.csv",
        "results": output_dir / "gene_level_results.csv.gz",
        "manifest": output_dir / "celltype_result_manifest.csv",
        "diagnostics": output_dir / "design_diagnostics.csv",
        "diagnostic_plots": output_dir / "edgeR_diagnostic_plots.pdf",
        "runtime_versions": output_dir / "runtime_versions.csv",
        "session_info": output_dir / "session_info.txt",
        "edger_log": output_dir / "edgeR.log",
    }
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
    manifest = pd.read_csv(paths["manifest"])
    diagnostics = pd.read_csv(paths["diagnostics"])
    versions = pd.read_csv(paths["runtime_versions"], dtype=str).fillna("")

    log_errors: list[str] = []
    log_contracts = [(paths["edger_log"], "return_code: 0")]
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

    profile_columns_ok = REQUIRED_PROFILE_COLUMNS.issubset(profiles.columns)
    profile_keys_ok = (
        profile_columns_ok
        and profiles["profile_id"].notna().all()
        and not profiles["profile_id"].duplicated().any()
        and not profiles[["biological_replicate_id", "cell_type"]].duplicated().any()
        and profiles[list(REQUIRED_PROFILE_COLUMNS)].notna().all().all()
        and profiles["n_cells"].ge(contract["min_cells_per_pseudobulk"]).all()
    )
    record(
        "profile_metadata_integrity",
        profile_keys_ok,
        f"Checked {len(profiles):,} unique replicate-by-cell-type profiles.",
    )

    primary = set(contract["primary_cell_types"])
    exploratory = set(contract["exploratory_cell_types"])
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
        (f"Observed {len(primary)} primary and {len(exploratory)} exploratory cell types."),
    )

    included = eligibility.loc[eligibility["aggregation_status"].eq("included"), "cell_type"]
    record(
        "eligibility_matches_contract",
        set(included) == approved and not included.duplicated().any(),
        f"Eligibility includes {len(included)} approved cell types.",
    )

    result_columns_ok = REQUIRED_RESULT_COLUMNS.issubset(results.columns)
    numeric_columns = [
        "log2_fc_per_10_years",
        "average_log_cpm",
        "ql_f_statistic",
        "p_value",
        "fdr_within_celltype",
        "fdr_global",
        "n_replicates",
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

    expected_tiers = {
        **{cell_type: "primary" for cell_type in primary},
        **{cell_type: "exploratory" for cell_type in exploratory},
    }
    result_tiers_ok = all(
        set(group["analysis_tier"]) == {expected_tiers[cell_type]}
        for cell_type, group in results.groupby("cell_type", sort=False)
    )
    model_formula_ok = (
        results.loc[results["model_scope"].eq("adjusted_primary"), "model_formula"]
        .eq(contract["primary_formula"])
        .all()
    )
    record(
        "result_analysis_contract",
        set(results["cell_type"]) == approved and result_tiers_ok and model_formula_ok,
        "Result tiers and adjusted formulas match the approved contract.",
    )

    diagnostic_cell_types = set(diagnostics["cell_type"])
    diagnostics_ok = (
        diagnostic_cell_types == approved
        and diagnostics["cell_type"].is_unique
        and diagnostics["status"].eq("completed").all()
        and diagnostics["design_rank"].eq(diagnostics["design_columns"]).all()
        and diagnostics["residual_df"].ge(contract["min_residual_df"]).all()
        and diagnostics["n_genes_tested"].sum() == len(results)
    )
    record(
        "model_diagnostics",
        diagnostics_ok,
        (
            "All approved cell types completed with full-rank designs and "
            "sufficient residual degrees of freedom."
        ),
    )

    fallback = diagnostics.loc[
        ~diagnostics["model_scope"].eq("adjusted_primary"), "cell_type"
    ].tolist()
    record(
        "adjusted_primary_models",
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
    )
    record(
        "runtime_versions",
        versions_ok,
        (
            f"Observed R {observed_versions.get('R')}, Bioconductor "
            f"{observed_versions.get('Bioconductor')}, and edgeR "
            f"{observed_versions.get('edgeR')}."
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
        "n_biological_replicates": int(profiles["biological_replicate_id"].nunique()),
        "n_genes_input": int(len(genes)),
        "n_cell_types_tested": int(results["cell_type"].nunique()),
        "n_primary_cell_types": int(len(primary)),
        "n_exploratory_cell_types": int(len(exploratory)),
        "n_gene_celltype_tests": int(len(results)),
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
