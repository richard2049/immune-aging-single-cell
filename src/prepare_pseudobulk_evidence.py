from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .utils import ensure_dir, load_config

PRIMARY_COLUMNS = {
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
SENSITIVITY_COLUMNS = {
    "gene_id",
    "cell_type",
    "analysis_tier",
    "scenario_id",
    "scenario_type",
    "excluded_batch",
    "model_formula",
    "log2_fc_per_10_years",
    "p_value",
    "n_replicates",
    "age_span_years",
}
DIAGNOSTIC_COLUMNS = {
    "cell_type",
    "analysis_tier",
    "scenario_id",
    "scenario_type",
    "model_formula",
    "n_replicates",
    "age_span_years",
    "design_columns",
    "design_rank",
    "residual_df",
    "n_genes_tested",
    "n_candidates_expected",
    "n_candidates_observed",
    "status",
    "reason",
}
PROFILE_COLUMNS = {
    "profile_id",
    "biological_replicate_id",
    "cell_type",
    "age",
    "sex",
    "batch",
    "analysis_tier",
    "n_cells",
    "library_size",
}
REFERENCE_DIAGNOSTIC_COLUMNS = {
    "cell_type",
    "analysis_tier",
    "n_replicates",
    "residual_df",
    "n_genes_tested",
}
SCENARIO_MANIFEST_COLUMNS = {
    "scenario_id",
    "scenario_type",
    "model_formula",
    "min_cells",
    "excluded_batch",
}
RUNTIME_VERSION_COLUMNS = {"component", "expected", "observed"}


def _require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{label} lacks required columns: {', '.join(missing)}")


def _write_csv(table: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_suffix(path.suffix + ".tmp")
    compression = "gzip" if path.suffix == ".gz" else None
    table.to_csv(temporary, index=False, compression=compression)
    temporary.replace(path)


def _write_json(payload: dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_suffix(path.suffix + ".tmp")

    def json_default(value: Any) -> Any:
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")

    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _direction(values: pd.Series | np.ndarray) -> np.ndarray:
    return np.sign(np.asarray(values, dtype=float)).astype(int)


def _safe_spearman(left: pd.Series, right: pd.Series) -> float:
    if left.nunique(dropna=True) < 2 or right.nunique(dropna=True) < 2:
        return float("nan")
    return float(left.corr(right, method="spearman"))


def _effect_for_scenario(table: pd.DataFrame, scenario_id: str) -> float:
    selected = table.loc[table["scenario_id"].eq(scenario_id), "log2_fc_per_10_years"]
    if len(selected) > 1:
        raise ValueError(f"Duplicate sensitivity result for {scenario_id}.")
    return float(selected.iloc[0]) if len(selected) else float("nan")


def _validate_sensitivity_contract(
    primary_results: pd.DataFrame,
    sensitivity_results: pd.DataFrame,
    scenario_diagnostics: pd.DataFrame,
    scenario_manifest: pd.DataFrame,
    global_fdr_threshold: float,
) -> None:
    candidates = primary_results.loc[
        primary_results["fdr_global"].lt(global_fdr_threshold),
        ["gene_id", "cell_type", "analysis_tier"],
    ].copy()
    if candidates[["gene_id", "cell_type"]].duplicated().any():
        raise ValueError("Primary candidate gene-by-cell-type keys are not unique.")

    candidate_tiers = candidates.rename(columns={"analysis_tier": "expected_tier"})
    observed = sensitivity_results.merge(
        candidate_tiers,
        on=["gene_id", "cell_type"],
        how="left",
        validate="many_to_one",
    )
    if observed["expected_tier"].isna().any():
        raise ValueError("Sensitivity output contains non-candidate gene keys.")
    if not observed["analysis_tier"].eq(observed["expected_tier"]).all():
        raise ValueError("Sensitivity result tiers do not match primary candidates.")

    if scenario_manifest["scenario_id"].duplicated().any():
        raise ValueError("Scenario manifest identifiers must be unique.")
    manifest_lookup = scenario_manifest.set_index("scenario_id")
    unknown_scenarios = set(scenario_diagnostics["scenario_id"]) - set(manifest_lookup.index)
    if unknown_scenarios:
        raise ValueError("Scenario diagnostics contain identifiers absent from the manifest.")

    candidate_counts = candidates.groupby("cell_type", observed=True).size()
    observed_counts = sensitivity_results.groupby(
        ["cell_type", "scenario_id"], observed=True
    ).size()
    diagnostic_celltypes = set(scenario_diagnostics["cell_type"])
    if diagnostic_celltypes != set(candidate_counts.index):
        raise ValueError("Scenario diagnostics do not cover every candidate cell type exactly.")

    for row in scenario_diagnostics.itertuples(index=False):
        manifest_row = manifest_lookup.loc[row.scenario_id]
        for field in ["scenario_type", "model_formula", "excluded_batch"]:
            if str(getattr(row, field)) != str(manifest_row[field]):
                raise ValueError(
                    f"Scenario diagnostic {row.scenario_id!r} disagrees with its manifest."
                )
        expected = int(candidate_counts.loc[row.cell_type])
        observed_count = int(observed_counts.get((row.cell_type, row.scenario_id), 0))
        if int(row.n_candidates_expected) != expected:
            raise ValueError("Scenario diagnostics contain an incorrect candidate denominator.")
        if int(row.n_candidates_observed) != observed_count:
            raise ValueError("Scenario diagnostics disagree with sensitivity result counts.")
        if row.status != "completed" and observed_count != 0:
            raise ValueError("Non-completed scenarios must not contain sensitivity results.")

    required_ids = {"high_cell_support", "drop_sex", "drop_batch"}
    if not required_ids.issubset(set(scenario_manifest["scenario_id"])):
        raise ValueError("Scenario manifest lacks a required approved sensitivity.")


def build_candidate_evidence(
    primary_results: pd.DataFrame,
    sensitivity_results: pd.DataFrame,
    scenario_diagnostics: pd.DataFrame,
    global_fdr_threshold: float,
) -> pd.DataFrame:
    candidates = primary_results.loc[primary_results["fdr_global"].lt(global_fdr_threshold)].copy()
    if candidates.empty:
        raise ValueError("No candidates pass the configured global FDR threshold.")
    if candidates[["gene_id", "cell_type"]].duplicated().any():
        raise ValueError("Primary candidate gene-by-cell-type keys are not unique.")

    candidate_keys = set(zip(candidates["gene_id"], candidates["cell_type"]))
    sensitivity_keys = set(zip(sensitivity_results["gene_id"], sensitivity_results["cell_type"]))
    unexpected = sensitivity_keys - candidate_keys
    if unexpected:
        raise ValueError("Sensitivity output contains non-candidate gene keys.")
    if sensitivity_results[["gene_id", "cell_type", "scenario_id"]].duplicated().any():
        raise ValueError("Sensitivity gene-by-cell-type-by-scenario keys are not unique.")

    completed = scenario_diagnostics.loc[scenario_diagnostics["status"].eq("completed")]
    completed_by_celltype = {
        cell_type: group.copy() for cell_type, group in completed.groupby("cell_type", sort=False)
    }
    sensitivity_by_key = {
        key: group.copy()
        for key, group in sensitivity_results.groupby(["gene_id", "cell_type"], sort=False)
    }

    records: list[dict[str, Any]] = []
    for row in candidates.itertuples(index=False):
        key = (row.gene_id, row.cell_type)
        observed = sensitivity_by_key.get(key, sensitivity_results.iloc[0:0])
        expected_diagnostics = completed_by_celltype.get(
            row.cell_type, scenario_diagnostics.iloc[0:0]
        )
        expected_scenarios = set(expected_diagnostics["scenario_id"])
        observed_scenarios = set(observed["scenario_id"])
        if not observed_scenarios.issubset(expected_scenarios):
            raise ValueError(f"Candidate {key!r} has results from a non-completed scenario.")

        primary_effect = float(row.log2_fc_per_10_years)
        effects = observed["log2_fc_per_10_years"].to_numpy(dtype=float)
        concordant = _direction(effects) == int(np.sign(primary_effect))
        leave_out = observed.loc[observed["scenario_type"].eq("leave_one_batch_out")]
        leave_effects = leave_out["log2_fc_per_10_years"].to_numpy(dtype=float)
        leave_concordant = _direction(leave_effects) == int(np.sign(primary_effect))

        records.append(
            {
                "gene_id": row.gene_id,
                "cell_type": row.cell_type,
                "analysis_tier": row.analysis_tier,
                "model_scope": row.model_scope,
                "model_formula": row.model_formula,
                "primary_log2_fc_per_10_years": primary_effect,
                "primary_average_log_cpm": float(row.average_log_cpm),
                "primary_ql_f_statistic": float(row.ql_f_statistic),
                "primary_p_value": float(row.p_value),
                "primary_fdr_within_celltype": float(row.fdr_within_celltype),
                "primary_fdr_global": float(row.fdr_global),
                "primary_n_replicates": int(row.n_replicates),
                "primary_age_span_years": float(row.age_span_years),
                "effect_direction": "higher_with_age" if primary_effect > 0 else "lower_with_age",
                "n_completed_scenarios": int(len(expected_scenarios)),
                "n_observed_scenarios": int(len(observed_scenarios)),
                "n_missing_after_expression_filter": int(
                    len(expected_scenarios - observed_scenarios)
                ),
                "direction_concordance_fraction": float(np.mean(concordant))
                if len(concordant)
                else float("nan"),
                "sensitivity_effect_min": float(np.min(effects)) if len(effects) else float("nan"),
                "sensitivity_effect_max": float(np.max(effects)) if len(effects) else float("nan"),
                "max_absolute_effect_change": float(np.max(np.abs(effects - primary_effect)))
                if len(effects)
                else float("nan"),
                "high_cell_support_effect": _effect_for_scenario(observed, "high_cell_support"),
                "drop_sex_effect": _effect_for_scenario(observed, "drop_sex"),
                "drop_batch_effect": _effect_for_scenario(observed, "drop_batch"),
                "n_leave_one_batch_out_completed": int(
                    expected_diagnostics["scenario_type"].eq("leave_one_batch_out").sum()
                ),
                "n_leave_one_batch_out_observed": int(len(leave_out)),
                "leave_one_batch_out_direction_concordance_fraction": float(
                    np.mean(leave_concordant)
                )
                if len(leave_concordant)
                else float("nan"),
                "leave_one_batch_out_effect_min": float(np.min(leave_effects))
                if len(leave_effects)
                else float("nan"),
                "leave_one_batch_out_effect_max": float(np.max(leave_effects))
                if len(leave_effects)
                else float("nan"),
                "requires_human_review": True,
                "automated_disposition": "not_assigned",
            }
        )
    evidence = pd.DataFrame.from_records(records)
    return evidence.sort_values(
        ["analysis_tier", "cell_type", "primary_fdr_global", "gene_id"],
        kind="stable",
    ).reset_index(drop=True)


def build_review_queue(evidence: pd.DataFrame, limit: int) -> pd.DataFrame:
    if limit <= 0:
        raise ValueError("Review candidates per direction and cell type must be positive.")
    selected: list[pd.DataFrame] = []
    for (_, _), group in evidence.groupby(["cell_type", "effect_direction"], sort=True):
        ranked = group.sort_values(
            ["primary_fdr_global", "max_absolute_effect_change", "gene_id"],
            ascending=[True, True, True],
            kind="stable",
        ).head(limit)
        ranked = ranked.copy()
        ranked["review_rank_within_direction"] = np.arange(1, len(ranked) + 1)
        ranked["selection_reason"] = (
            "top_global_fdr_within_celltype_and_direction_for_manual_review"
        )
        selected.append(ranked)
    if not selected:
        return evidence.iloc[0:0].copy()
    return pd.concat(selected, ignore_index=True).sort_values(
        ["analysis_tier", "cell_type", "effect_direction", "review_rank_within_direction"],
        kind="stable",
    )


def _load_diagnostic_review(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Diagnostic review record not found: {path}")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Diagnostic review record must be a YAML mapping.")
    return payload


def build_celltype_summary(
    profiles: pd.DataFrame,
    reference_diagnostics: pd.DataFrame,
    scenario_diagnostics: pd.DataFrame,
    evidence: pd.DataFrame,
    review_queue: pd.DataFrame,
    diagnostic_review: dict[str, Any],
) -> pd.DataFrame:
    page_records = diagnostic_review.get("pages", [])
    page_status = {
        str(record["cell_type"]): str(record["status"])
        for record in page_records
        if isinstance(record, dict) and "cell_type" in record and "status" in record
    }
    records: list[dict[str, Any]] = []
    for row in reference_diagnostics.itertuples(index=False):
        cell_type = row.cell_type
        profile_group = profiles.loc[profiles["cell_type"].eq(cell_type)].copy()
        scenario_group = scenario_diagnostics.loc[scenario_diagnostics["cell_type"].eq(cell_type)]
        candidate_group = evidence.loc[evidence["cell_type"].eq(cell_type)]
        batch_counts = profile_group.groupby("batch", observed=True).size()
        records.append(
            {
                "cell_type": cell_type,
                "analysis_tier": row.analysis_tier,
                "reference_n_replicates": int(row.n_replicates),
                "reference_residual_df": int(row.residual_df),
                "reference_n_genes_tested": int(row.n_genes_tested),
                "n_global_fdr_candidates": int(len(candidate_group)),
                "n_review_queue_candidates": int(review_queue["cell_type"].eq(cell_type).sum()),
                "n_profiles": int(len(profile_group)),
                "n_batches": int(profile_group["batch"].nunique()),
                "minimum_profiles_per_batch": int(batch_counts.min()),
                "n_sex_levels": int(profile_group["sex"].nunique()),
                "minimum_cells_per_profile": int(profile_group["n_cells"].min()),
                "median_cells_per_profile": float(profile_group["n_cells"].median()),
                "minimum_library_size": int(profile_group["library_size"].min()),
                "median_library_size": float(profile_group["library_size"].median()),
                "spearman_age_vs_n_cells": _safe_spearman(
                    profile_group["age"], profile_group["n_cells"]
                ),
                "spearman_age_vs_library_size": _safe_spearman(
                    profile_group["age"], profile_group["library_size"]
                ),
                "n_scenarios_completed": int(scenario_group["status"].eq("completed").sum()),
                "n_scenarios_not_estimable": int(
                    scenario_group["status"].eq("not_estimable").sum()
                ),
                "n_scenarios_failed": int(scenario_group["status"].eq("model_failed").sum()),
                "minimum_candidate_direction_concordance": float(
                    candidate_group["direction_concordance_fraction"].min()
                )
                if len(candidate_group)
                else float("nan"),
                "median_candidate_direction_concordance": float(
                    candidate_group["direction_concordance_fraction"].median()
                )
                if len(candidate_group)
                else float("nan"),
                "diagnostic_plot_review_status": page_status.get(cell_type, "missing"),
                "requires_human_interpretation": True,
            }
        )
    return pd.DataFrame.from_records(records).sort_values(
        ["analysis_tier", "cell_type"], kind="stable"
    )


def _render_report(
    summary: dict[str, Any],
    celltypes: pd.DataFrame,
    scenario_diagnostics: pd.DataFrame,
    diagnostic_review: dict[str, Any],
) -> str:
    execution_labels = {
        "completed_log_and_validated_outputs": "complete log and validated outputs",
        "validated_outputs_after_incomplete_wrapper_log": (
            "validated outputs; wrapper log incomplete after Docker client hang"
        ),
    }
    diagnostic_labels = {
        "completed_with_followup_flags": "complete; follow-up flags recorded",
        "reviewed_no_blocking_issue": "reviewed; no blocking issue",
        "reviewed_followup_required": "reviewed; follow-up required",
    }
    status_counts = (
        scenario_diagnostics.groupby(["scenario_type", "status"], dropna=False)
        .size()
        .reset_index(name="models")
    )
    lines = [
        "# Pseudobulk Robustness Audit And Evidence Preparation",
        "",
        "## Scope",
        "",
        "This stage evaluates technical and statistical stability of the frozen",
        "donor-aware pseudobulk results. It does not approve genes, pathways,",
        "biological interpretations, or public claims.",
        "",
        "## Audit Status",
        "",
        f"- Automated contract: {'passed' if summary['passed'] else 'failed'}",
        (
            f"- Global-FDR candidate rows preserved: {summary['n_candidates']:,} "
            f"({summary['n_primary_candidates']:,} primary; "
            f"{summary['n_exploratory_candidates']:,} exploratory)"
        ),
        f"- Compact manual-review queue: {summary['n_review_queue']:,}",
        f"- Completed sensitivity models: {summary['n_completed_models']:,}",
        "- Sensitivity execution record: "
        + execution_labels.get(
            summary["execution_record_status"], summary["execution_record_status"]
        ),
        "- Diagnostic plot review: "
        + diagnostic_labels.get(
            str(diagnostic_review.get("overall_status", "missing")),
            str(diagnostic_review.get("overall_status", "missing")),
        ),
        "- Biological interpretation: not performed",
        "",
        "## Sensitivity Model Inventory",
        "",
        "| Scenario type | Status | Models |",
        "| --- | --- | ---: |",
    ]
    for row in status_counts.itertuples(index=False):
        scenario_type = str(row.scenario_type).replace("_", " ")
        status = str(row.status).replace("_", " ")
        lines.append(f"| {scenario_type} | {status} | {row.models} |")
    non_estimable = scenario_diagnostics.loc[scenario_diagnostics["status"].eq("not_estimable")]
    if not non_estimable.empty:
        lines.extend(
            [
                "",
                "## Non-Estimable Scenarios",
                "",
                "| Cell type | Tier | Scenario | Replicates | Residual df | Reason |",
                "| --- | --- | --- | ---: | ---: | --- |",
            ]
        )
        for row in non_estimable.itertuples(index=False):
            scenario = str(row.scenario_id).replace("__", ": ").replace("_", " ")
            reason = str(row.reason).replace("_", " ")
            residual_df = "NA" if not np.isfinite(row.residual_df) else str(int(row.residual_df))
            lines.append(
                f"| {row.cell_type} | {row.analysis_tier} | {scenario} | "
                f"{row.n_replicates} | {residual_df} | {reason} |"
            )
    lines.extend(
        [
            "",
            "## Cell-Type Audit Summary",
            "",
            "| Cell type | Tier | Replicates | Residual df | Global-FDR candidates | Completed scenarios | Median direction concordance | MD plot review |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in celltypes.itertuples(index=False):
        concordance = (
            "NA"
            if not np.isfinite(row.median_candidate_direction_concordance)
            else f"{row.median_candidate_direction_concordance:.3f}"
        )
        lines.append(
            f"| {row.cell_type} | {row.analysis_tier} | "
            f"{row.reference_n_replicates} | {row.reference_residual_df} | "
            f"{row.n_global_fdr_candidates} | {row.n_scenarios_completed} | "
            f"{concordance} | "
            f"{diagnostic_labels.get(row.diagnostic_plot_review_status, row.diagnostic_plot_review_status)} |"
        )
    lines.extend(
        [
            "",
            "## Runtime Provenance",
            "",
            "| Component | Expected | Observed |",
            "| --- | --- | --- |",
        ]
    )
    for component, versions in summary["runtime_versions"].items():
        lines.append(
            f"| {component} | {versions['expected'] or 'not pinned'} | {versions['observed']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "Sensitivity models reuse the same cohort and are not independent",
            "replications. Direction concordance and effect ranges are review",
            "diagnostics, not automated evidence of biological validity. Omission",
            "of sex or batch deliberately weakens the adjustment set and is used",
            "only to expose adjustment dependence. Leave-one-batch-out models can",
            "identify sensitivity to a processing batch but cannot remove all",
            "residual confounding.",
            "",
            "The exploratory NK-cell population remains separate because its",
            "reference model used only 18 qualifying replicates and eight residual",
            "degrees of freedom. Gene annotation, pathway enrichment, candidate",
            "disposition, external replication, and public wording require a later",
            "D-class review.",
            "",
        ]
    )
    if summary["warnings"]:
        lines.extend(["## Operational Notes", ""])
        lines.extend(f"- {warning}" for warning in summary["warnings"])
        lines.append("")
    return "\n".join(lines)


def prepare(
    config_path: str | Path,
    pseudobulk_dir: str | Path,
    report_path: str | Path,
) -> dict[str, Any]:
    config_path = Path(config_path)
    pseudobulk_dir = Path(pseudobulk_dir)
    report_path = Path(report_path)
    config = load_config(config_path)
    robustness = config["pseudobulk_robustness"]
    threshold = float(robustness["candidate_global_fdr_threshold"])
    queue_limit = int(robustness["review_candidates_per_direction_per_celltype"])
    if not 0 < threshold < 1:
        raise ValueError("Candidate global FDR threshold must be between zero and one.")

    robustness_dir = pseudobulk_dir / "robustness"
    primary_results = pd.read_csv(pseudobulk_dir / "gene_level_results.csv.gz")
    profiles = pd.read_csv(pseudobulk_dir / "pseudobulk_profiles.csv")
    reference_diagnostics = pd.read_csv(pseudobulk_dir / "design_diagnostics.csv")
    sensitivity_results = pd.read_csv(robustness_dir / "sensitivity_gene_results.csv.gz")
    scenario_diagnostics = pd.read_csv(robustness_dir / "scenario_diagnostics.csv").fillna(
        {"reason": "", "excluded_batch": ""}
    )
    scenario_manifest = pd.read_csv(robustness_dir / "scenario_manifest.csv").fillna(
        {"excluded_batch": ""}
    )
    technical_validation = json.loads(
        (pseudobulk_dir / "technical_validation.json").read_text(encoding="utf-8")
    )
    if not bool(technical_validation.get("passed")):
        raise ValueError("Primary pseudobulk technical validation has not passed.")
    runtime_versions = pd.read_csv(pseudobulk_dir / "runtime_versions.csv").fillna("")
    robustness_log = robustness_dir / "edgeR_robustness.log"
    if not robustness_log.is_file():
        raise FileNotFoundError(f"Robustness execution log not found: {robustness_log}")
    robustness_log_text = robustness_log.read_text(encoding="utf-8")
    diagnostic_review_path = Path(robustness["diagnostic_review_path"])
    diagnostic_review = _load_diagnostic_review(diagnostic_review_path)

    _require_columns(primary_results, PRIMARY_COLUMNS, "Primary results")
    _require_columns(sensitivity_results, SENSITIVITY_COLUMNS, "Sensitivity results")
    _require_columns(scenario_diagnostics, DIAGNOSTIC_COLUMNS, "Scenario diagnostics")
    _require_columns(profiles, PROFILE_COLUMNS, "Pseudobulk profiles")
    _require_columns(
        reference_diagnostics,
        REFERENCE_DIAGNOSTIC_COLUMNS,
        "Reference diagnostics",
    )
    _require_columns(scenario_manifest, SCENARIO_MANIFEST_COLUMNS, "Scenario manifest")
    _require_columns(runtime_versions, RUNTIME_VERSION_COLUMNS, "Runtime versions")
    if scenario_diagnostics[["cell_type", "scenario_id"]].duplicated().any():
        raise ValueError("Scenario diagnostics contain duplicate cell-type keys.")
    if reference_diagnostics["cell_type"].duplicated().any():
        raise ValueError("Reference diagnostics contain duplicate cell types.")
    _validate_sensitivity_contract(
        primary_results,
        sensitivity_results,
        scenario_diagnostics,
        scenario_manifest,
        threshold,
    )

    evidence = build_candidate_evidence(
        primary_results,
        sensitivity_results,
        scenario_diagnostics,
        threshold,
    )
    review_queue = build_review_queue(evidence, queue_limit)
    celltype_summary = build_celltype_summary(
        profiles,
        reference_diagnostics,
        scenario_diagnostics,
        evidence,
        review_queue,
        diagnostic_review,
    )

    expected_candidates = primary_results.loc[
        primary_results["fdr_global"].lt(threshold), ["gene_id", "cell_type"]
    ]
    observed_candidates = evidence[["gene_id", "cell_type"]]
    review_pages = diagnostic_review.get("pages", [])
    reviewed_celltypes = {
        str(record.get("cell_type")) for record in review_pages if isinstance(record, dict)
    }
    expected_celltypes = set(reference_diagnostics["cell_type"])
    valid_review_statuses = {"reviewed_no_blocking_issue", "reviewed_followup_required"}
    page_statuses = {
        str(record.get("status")) for record in review_pages if isinstance(record, dict)
    }
    log_completed = (
        "return_code: 0" in robustness_log_text
        and "[pseudobulk_robustness] completed" in robustness_log_text
    )
    execution_record_status = (
        "completed_log_and_validated_outputs"
        if log_completed
        else "validated_outputs_after_incomplete_wrapper_log"
    )
    warnings = []
    if not log_completed:
        warnings.append(
            "The Docker container wrote the complete sensitivity outputs, but "
            "the Windows Docker client did not close its wrapper log. The two "
            "identified host processes were stopped after the container exited; "
            "acceptance relies on independent schema, key, count, scenario, and "
            "model-status checks."
        )

    checks = [
        {
            "name": "technical_validation_passed",
            "passed": True,
        },
        {
            "name": "candidate_set_preserved",
            "passed": set(map(tuple, expected_candidates.to_numpy()))
            == set(map(tuple, observed_candidates.to_numpy())),
        },
        {
            "name": "no_automated_biological_dispositions",
            "passed": evidence["automated_disposition"].eq("not_assigned").all()
            and evidence["requires_human_review"].all(),
        },
        {
            "name": "review_queue_is_bounded",
            "passed": bool(
                review_queue.groupby(["cell_type", "effect_direction"]).size().le(queue_limit).all()
            ),
        },
        {
            "name": "diagnostic_plot_pages_reviewed",
            "passed": reviewed_celltypes == expected_celltypes
            and page_statuses.issubset(valid_review_statuses)
            and len(review_pages) == len(expected_celltypes),
        },
        {
            "name": "scenario_manifest_nonempty",
            "passed": not scenario_manifest.empty,
        },
        {
            "name": "sensitivity_output_contract_complete",
            "passed": True,
        },
        {
            "name": "no_sensitivity_model_failures",
            "passed": not scenario_diagnostics["status"].eq("model_failed").any(),
        },
    ]
    passed = all(check["passed"] for check in checks)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path),
        "pseudobulk_dir": str(pseudobulk_dir),
        "risk_classification": "C",
        "passed": passed,
        "biological_interpretation_status": "not_performed",
        "requires_d_class_human_review": True,
        "candidate_global_fdr_threshold": threshold,
        "n_candidates": int(len(evidence)),
        "n_primary_candidates": int(evidence["analysis_tier"].eq("primary").sum()),
        "n_exploratory_candidates": int(evidence["analysis_tier"].eq("exploratory").sum()),
        "n_review_queue": int(len(review_queue)),
        "n_completed_models": int(scenario_diagnostics["status"].eq("completed").sum()),
        "n_not_estimable_models": int(scenario_diagnostics["status"].eq("not_estimable").sum()),
        "n_failed_models": int(scenario_diagnostics["status"].eq("model_failed").sum()),
        "execution_record_status": execution_record_status,
        "warnings": warnings,
        "runtime_versions": {
            str(row.component): {
                "expected": str(row.expected),
                "observed": str(row.observed),
            }
            for row in runtime_versions.itertuples(index=False)
        },
        "checks": checks,
    }

    _write_csv(evidence, robustness_dir / "candidate_evidence.csv.gz")
    _write_csv(review_queue, robustness_dir / "review_queue.csv")
    _write_csv(
        celltype_summary,
        robustness_dir / "celltype_robustness_summary.csv",
    )
    _write_json(summary, robustness_dir / "audit_summary.json")
    ensure_dir(report_path.parent)
    report_path.write_text(
        _render_report(
            summary,
            celltype_summary,
            scenario_diagnostics,
            diagnostic_review,
        ),
        encoding="utf-8",
    )
    if not passed:
        failed = [check["name"] for check in checks if not check["passed"]]
        raise ValueError("Pseudobulk robustness audit failed: " + ", ".join(failed))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare pseudobulk robustness evidence without assigning biological "
            "interpretations or claims."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--pseudobulk-dir", required=True)
    parser.add_argument("--report-out", required=True)
    args = parser.parse_args()
    summary = prepare(args.config, args.pseudobulk_dir, args.report_out)
    print(
        "[prepare_pseudobulk_evidence] passed: "
        f"{summary['n_candidates']:,} candidates, "
        f"{summary['n_review_queue']:,} queued for human review"
    )


if __name__ == "__main__":
    main()
