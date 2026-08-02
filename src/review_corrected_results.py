from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.utils import ensure_dir, load_config


ASSOCIATION_GROUP = "biological_replicate_id"
ANALYSIS_SPECS = {
    "composition": {
        "key": ["cell_type"],
        "effect": "slope_per_10y",
        "ci_low": "slope_per_10y_ci_low",
        "ci_high": "slope_per_10y_ci_high",
        "association": "spearman_rho",
        "association_ci_low": "spearman_rho_ci_low",
        "association_ci_high": "spearman_rho_ci_high",
        "fdr": "spearman_fdr",
    },
    "signature": {
        "key": ["cell_type", "signature"],
        "effect": "effect_per_10y",
        "ci_low": "effect_per_10y_ci_low",
        "ci_high": "effect_per_10y_ci_high",
        "association": "spearman_rho",
        "association_ci_low": "spearman_rho_ci_low",
        "association_ci_high": "spearman_rho_ci_high",
        "fdr": "fdr",
    },
}
ALLOWED_DISPOSITIONS = {"retain", "exploratory", "exclude", "defer"}
ALLOWED_PRESENTATION_TIERS = {
    "core",
    "secondary",
    "exploratory",
    "excluded",
    "technical_only",
}


def _require_columns(
    frame: pd.DataFrame,
    columns: list[str],
    source: Path,
) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def _as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def _finite_float(value: Any) -> float:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(number) if pd.notna(number) else np.nan


def _clean_text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def _effect_direction(value: Any) -> int:
    number = _finite_float(value)
    if not np.isfinite(number) or number == 0:
        return 0
    return 1 if number > 0 else -1


def _ci_excludes_zero(low: Any, high: Any) -> bool:
    low_value = _finite_float(low)
    high_value = _finite_float(high)
    return bool(
        np.isfinite(low_value)
        and np.isfinite(high_value)
        and (low_value > 0 or high_value < 0)
    )


def _signal_id(
    analysis: str,
    cell_type: Any = "",
    signature: Any = "",
) -> str:
    parts = [analysis, str(cell_type).strip()]
    if str(signature).strip():
        parts.append(str(signature).strip())
    return "::".join(parts)


def _load_comparison(path: Path) -> dict[tuple[str, str, str], dict]:
    comparison = pd.read_csv(path)
    _require_columns(
        comparison,
        ["analysis", "cell_type", "comparison_status"],
        path,
    )
    records: dict[tuple[str, str, str], dict] = {}
    for row in comparison.to_dict("records"):
        key = (
            _clean_text(row["analysis"]),
            _clean_text(row["cell_type"]),
            _clean_text(row.get("signature", "")),
        )
        records[key] = row
    return records


def _composition_support(path: Path) -> dict[str, dict[str, float]]:
    frame = pd.read_csv(path)
    _require_columns(
        frame,
        [
            ASSOCIATION_GROUP,
            "age",
            "cell_type",
            "n_cells",
            "total_cells",
        ],
        path,
    )
    support: dict[str, dict[str, float]] = {}
    for cell_type, group in frame.groupby("cell_type", observed=True):
        eligible = group.loc[group["n_cells"].gt(0)]
        support[str(cell_type)] = {
            "age_min": _finite_float(eligible["age"].min()),
            "age_max": _finite_float(eligible["age"].max()),
            "age_span": _finite_float(
                eligible["age"].max() - eligible["age"].min()
            ),
            "total_cells": int(eligible["n_cells"].sum()),
            "median_cells_per_replicate": _finite_float(
                eligible["n_cells"].median()
            ),
        }
    return support


def _load_sensitivity(
    manifest_path: Path,
) -> tuple[pd.DataFrame, dict[str, dict[str, list[dict]]]]:
    manifest = pd.read_csv(manifest_path)
    _require_columns(
        manifest,
        [
            "scenario",
            "composition_trends_table",
            "signature_assoc_table",
        ],
        manifest_path,
    )
    evidence: dict[str, dict[str, list[dict]]] = {
        "composition": {},
        "signature": {},
    }
    for scenario in manifest.to_dict("records"):
        for analysis, table_column in (
            ("composition", "composition_trends_table"),
            ("signature", "signature_assoc_table"),
        ):
            table_path = Path(str(scenario[table_column]))
            table = pd.read_csv(table_path)
            spec = ANALYSIS_SPECS[analysis]
            required = (
                spec["key"]
                + [
                    spec["effect"],
                    spec["fdr"],
                    "n_donors",
                    "grouping_id_column",
                ]
            )
            _require_columns(table, required, table_path)
            for row in table.to_dict("records"):
                signal = _signal_id(
                    analysis,
                    row["cell_type"],
                    row.get("signature", ""),
                )
                evidence[analysis].setdefault(signal, []).append(
                    {
                        "scenario": str(scenario["scenario"]),
                        "effect": _finite_float(row[spec["effect"]]),
                        "fdr": _finite_float(row[spec["fdr"]]),
                        "n_replicates": int(row["n_donors"]),
                        "grouping_id_column": str(
                            row["grouping_id_column"]
                        ),
                        "source_path": str(table_path),
                    }
                )
    return manifest, evidence


def _sensitivity_summary(
    records: list[dict],
    expected_scenarios: int,
    corrected_direction: int,
    fdr_threshold: float,
) -> dict[str, Any]:
    directions = [_effect_direction(item["effect"]) for item in records]
    same_direction = [
        direction == corrected_direction
        for direction in directions
        if direction != 0 and corrected_direction != 0
    ]
    valid_fdr = [
        item["fdr"]
        for item in records
        if np.isfinite(item["fdr"])
    ]
    effects = [
        item["effect"]
        for item in records
        if np.isfinite(item["effect"])
    ]
    grouping_ok = all(
        item["grouping_id_column"] == ASSOCIATION_GROUP
        for item in records
    )
    tested = len(records)
    direction_fraction = (
        float(np.mean(same_direction)) if same_direction else 0.0
    )
    fdr_count = sum(value < fdr_threshold for value in valid_fdr)
    return {
        "sensitivity_expected_scenarios": expected_scenarios,
        "sensitivity_tested_scenarios": tested,
        "sensitivity_direction_consistency_fraction": direction_fraction,
        "sensitivity_fdr_supported_count": fdr_count,
        "sensitivity_fdr_supported_fraction": (
            fdr_count / expected_scenarios if expected_scenarios else np.nan
        ),
        "sensitivity_effect_min": min(effects) if effects else np.nan,
        "sensitivity_effect_max": max(effects) if effects else np.nan,
        "sensitivity_min_replicates": min(
            (item["n_replicates"] for item in records),
            default=0,
        ),
        "sensitivity_grouping_valid": grouping_ok,
        "sensitivity_complete_pass": tested == expected_scenarios,
        "direction_consistency_pass": (
            tested == expected_scenarios and all(same_direction)
        ),
        "all_sensitivity_fdr_pass": (
            tested == expected_scenarios
            and len(valid_fdr) == expected_scenarios
            and fdr_count == expected_scenarios
        ),
    }


def classify_association(
    row: dict[str, Any],
    criteria: dict[str, Any],
) -> str:
    core_pass = (
        row["corrected_fdr_pass"]
        and row["min_support_pass"]
        and (
            row["effect_ci_excludes_zero"]
            or not criteria["require_effect_ci_excludes_zero"]
        )
        and row["grouping_valid"]
    )
    if not core_pass:
        return "not_prioritized"

    sensitivity_pass = (
        (
            row["sensitivity_complete_pass"]
            or not criteria["require_all_sensitivity_scenarios"]
        )
        and (
            row["direction_consistency_pass"]
            or not criteria["require_direction_consistency"]
        )
        and (
            row["all_sensitivity_fdr_pass"]
            or not criteria["require_fdr_support_all_scenarios"]
        )
        and row["sensitivity_grouping_valid"]
    )
    if sensitivity_pass:
        return "candidate_for_human_review"
    return "exploratory_sensitivity_limited"


def _association_rows(
    analysis: str,
    source_path: Path,
    support: dict[str, dict[str, float]],
    sensitivity: dict[str, list[dict]],
    comparison: dict[tuple[str, str, str], dict],
    criteria: dict[str, Any],
    expected_scenarios: int,
) -> list[dict[str, Any]]:
    frame = pd.read_csv(source_path)
    spec = ANALYSIS_SPECS[analysis]
    required = (
        spec["key"]
        + [
            "n_donors",
            spec["effect"],
            spec["ci_low"],
            spec["ci_high"],
            spec["association"],
            spec["association_ci_low"],
            spec["association_ci_high"],
            spec["fdr"],
            "r_squared",
            "direction",
            "grouping_id_column",
        ]
    )
    _require_columns(frame, required, source_path)

    evidence_rows: list[dict[str, Any]] = []
    for source_row in frame.to_dict("records"):
        cell_type = str(source_row["cell_type"])
        signature = _clean_text(source_row.get("signature", ""))
        signal = _signal_id(analysis, cell_type, signature)
        effect = _finite_float(source_row[spec["effect"]])
        fdr = _finite_float(source_row[spec["fdr"]])
        support_row = support.get(cell_type, {})
        comparison_row = comparison.get(
            (analysis, cell_type, signature),
            {},
        )
        row = {
            "analysis": analysis,
            "signal_id": signal,
            "cell_type": cell_type,
            "signature": signature,
            "n_replicates": int(source_row["n_donors"]),
            "age_min": _finite_float(
                source_row.get("age_min", support_row.get("age_min"))
            ),
            "age_max": _finite_float(
                source_row.get("age_max", support_row.get("age_max"))
            ),
            "age_span": _finite_float(
                source_row.get("age_span", support_row.get("age_span"))
            ),
            "total_cells": support_row.get("total_cells", np.nan),
            "median_cells_per_replicate": support_row.get(
                "median_cells_per_replicate",
                np.nan,
            ),
            "effect_per_10y": effect,
            "effect_ci_low": _finite_float(source_row[spec["ci_low"]]),
            "effect_ci_high": _finite_float(source_row[spec["ci_high"]]),
            "association": _finite_float(
                source_row[spec["association"]]
            ),
            "association_ci_low": _finite_float(
                source_row[spec["association_ci_low"]]
            ),
            "association_ci_high": _finite_float(
                source_row[spec["association_ci_high"]]
            ),
            "primary_fdr": fdr,
            "r_squared": _finite_float(source_row["r_squared"]),
            "direction": str(source_row["direction"]),
            "grouping_id_column": str(
                source_row["grouping_id_column"]
            ),
            "grouping_valid": (
                str(source_row["grouping_id_column"])
                == ASSOCIATION_GROUP
            ),
            "effect_ci_excludes_zero": _ci_excludes_zero(
                source_row[spec["ci_low"]],
                source_row[spec["ci_high"]],
            ),
            "corrected_fdr_pass": (
                np.isfinite(fdr) and fdr < criteria["fdr_threshold"]
            ),
            "min_support_pass": (
                int(source_row["n_donors"])
                >= criteria["min_biological_replicates"]
            ),
            "comparison_status": str(
                comparison_row.get("comparison_status", "")
            ),
            "provisional_effect": _finite_float(
                comparison_row.get("provisional_effect")
            ),
            "provisional_fdr": _finite_float(
                comparison_row.get("provisional_fdr")
            ),
            "source_path": str(source_path),
            "source_row_key": signal,
        }
        row.update(
            _sensitivity_summary(
                sensitivity.get(signal, []),
                expected_scenarios,
                _effect_direction(effect),
                criteria["fdr_threshold"],
            )
        )
        row["automated_screen_status"] = classify_association(
            row,
            criteria,
        )
        row["human_disposition"] = "pending"
        row["public_claim_approved"] = False
        row["primary_limitation"] = (
            "Cross-sectional association; residual confounding and "
            "annotation uncertainty remain."
            if analysis == "composition"
            else "Configured gene-set proxy; it does not establish pathway "
            "activation, causality, or external validity."
        )
        evidence_rows.append(row)
    return evidence_rows


def _prediction_row(
    metrics_path: Path,
    comparison_path: Path,
    validation_path: Path,
) -> dict[str, Any]:
    metrics = pd.read_csv(metrics_path)
    _require_columns(
        metrics,
        [
            "model",
            "mae",
            "mae_ci_low",
            "mae_ci_high",
            "r2",
            "r2_ci_low",
            "r2_ci_high",
            "baseline_mae",
            "delta_mae_vs_baseline",
            "n_donors",
            "n_groups",
            "n_splits",
            "is_best",
            "grouping_id_column",
        ],
        metrics_path,
    )
    best = metrics.loc[metrics["is_best"].map(_as_bool)]
    if len(best) != 1:
        raise ValueError(
            f"{metrics_path} must contain exactly one is_best row; "
            f"found {len(best)}"
        )
    selected = best.iloc[0]

    prediction_comparison = pd.read_csv(comparison_path)
    comparison_record = {}
    if not prediction_comparison.empty:
        candidates = prediction_comparison.loc[
            prediction_comparison["model"].eq(selected["model"])
        ]
        if not candidates.empty:
            comparison_record = candidates.iloc[0].to_dict()

    with validation_path.open("r", encoding="utf-8") as handle:
        validation = json.load(handle)
    validation_passed = bool(validation.get("passed", False))
    grouping_valid = selected["grouping_id_column"] == ASSOCIATION_GROUP
    signal_detected = (
        _finite_float(selected["delta_mae_vs_baseline"]) > 0
        and _finite_float(selected["r2_ci_low"]) > 0
    )
    screen_status = (
        "candidate_internal_performance_only"
        if validation_passed and grouping_valid and signal_detected
        else "not_prioritized"
    )
    return {
        "analysis": "prediction",
        "signal_id": f"prediction::{selected['model']}",
        "cell_type": "",
        "signature": "",
        "n_replicates": int(selected["n_donors"]),
        "age_min": np.nan,
        "age_max": np.nan,
        "age_span": np.nan,
        "total_cells": np.nan,
        "median_cells_per_replicate": np.nan,
        "effect_per_10y": np.nan,
        "effect_ci_low": np.nan,
        "effect_ci_high": np.nan,
        "association": _finite_float(selected.get("spearman_rho")),
        "association_ci_low": np.nan,
        "association_ci_high": np.nan,
        "primary_fdr": np.nan,
        "r_squared": _finite_float(selected["r2"]),
        "direction": "",
        "grouping_id_column": str(selected["grouping_id_column"]),
        "grouping_valid": grouping_valid,
        "effect_ci_excludes_zero": False,
        "corrected_fdr_pass": False,
        "min_support_pass": True,
        "comparison_status": str(
            comparison_record.get("comparison_status", "")
        ),
        "provisional_effect": np.nan,
        "provisional_fdr": np.nan,
        "source_path": str(metrics_path),
        "source_row_key": f"model={selected['model']}",
        "sensitivity_expected_scenarios": np.nan,
        "sensitivity_tested_scenarios": np.nan,
        "sensitivity_direction_consistency_fraction": np.nan,
        "sensitivity_fdr_supported_count": np.nan,
        "sensitivity_fdr_supported_fraction": np.nan,
        "sensitivity_effect_min": np.nan,
        "sensitivity_effect_max": np.nan,
        "sensitivity_min_replicates": np.nan,
        "sensitivity_grouping_valid": True,
        "sensitivity_complete_pass": False,
        "direction_consistency_pass": False,
        "all_sensitivity_fdr_pass": False,
        "automated_screen_status": screen_status,
        "human_disposition": "pending",
        "public_claim_approved": False,
        "primary_limitation": (
            "Internal grouped cross-validation only; no external validation "
            "and no evidence for clinical biomarker or aging-clock use."
        ),
        "prediction_model": str(selected["model"]),
        "mae_years": _finite_float(selected["mae"]),
        "mae_ci_low": _finite_float(selected["mae_ci_low"]),
        "mae_ci_high": _finite_float(selected["mae_ci_high"]),
        "r2_ci_low": _finite_float(selected["r2_ci_low"]),
        "r2_ci_high": _finite_float(selected["r2_ci_high"]),
        "baseline_mae_years": _finite_float(selected["baseline_mae"]),
        "delta_mae_vs_baseline": _finite_float(
            selected["delta_mae_vs_baseline"]
        ),
        "n_groups": int(selected["n_groups"]),
        "n_splits": int(selected["n_splits"]),
        "validation_passed": validation_passed,
    }


def apply_human_dispositions(
    evidence: pd.DataFrame,
    dispositions_path: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    record = load_config(dispositions_path)
    metadata = record.get("review_metadata", {})
    decisions = record.get("decisions", [])
    asset_decisions = record.get("draft_asset_decisions", {})
    final_approval = record.get("final_approval", {})
    if not isinstance(metadata, dict):
        raise TypeError(
            f"review_metadata must be a mapping: {dispositions_path}"
        )
    if not isinstance(decisions, list):
        raise TypeError(f"decisions must be a list: {dispositions_path}")
    if not isinstance(asset_decisions, dict):
        raise TypeError(
            f"draft_asset_decisions must be a mapping: {dispositions_path}"
        )
    if not isinstance(final_approval, dict):
        raise TypeError(
            f"final_approval must be a mapping: {dispositions_path}"
        )
    required_metadata = {
        "status",
        "reviewer",
        "review_date",
        "public_scope",
        "public_claims_approved",
        "asset_replacements_approved",
    }
    missing_metadata = sorted(required_metadata.difference(metadata))
    if missing_metadata:
        raise ValueError(
            f"{dispositions_path} is missing review metadata: "
            f"{missing_metadata}"
        )
    public_approved = _as_bool(metadata["public_claims_approved"])
    assets_approved = _as_bool(metadata["asset_replacements_approved"])
    final_status = str(final_approval.get("status", ""))
    if (public_approved or assets_approved) and final_status != "approved":
        raise ValueError(
            "Public claims or assets require final_approval.status=approved."
        )
    if public_approved and not final_approval.get("readme_wording"):
        raise ValueError("Public claims require approved README wording.")

    reviewed = evidence.copy()
    reviewed["presentation_tier"] = ""
    reviewed["disposition_rationale"] = ""
    reviewed["decision_reviewer"] = ""
    reviewed["decision_date"] = ""
    seen: set[str] = set()
    for decision in decisions:
        if not isinstance(decision, dict):
            raise TypeError(
                f"Each decision must be a mapping: {dispositions_path}"
            )
        required = {
            "signal_id",
            "disposition",
            "presentation_tier",
            "rationale",
        }
        missing = sorted(required.difference(decision))
        if missing:
            raise ValueError(
                f"A decision in {dispositions_path} is missing: {missing}"
            )
        signal_id = str(decision["signal_id"])
        if signal_id in seen:
            raise ValueError(f"Duplicate disposition for {signal_id}")
        seen.add(signal_id)
        disposition = str(decision["disposition"])
        tier = str(decision["presentation_tier"])
        if disposition not in ALLOWED_DISPOSITIONS:
            raise ValueError(
                f"Unsupported disposition for {signal_id}: {disposition}"
            )
        if tier not in ALLOWED_PRESENTATION_TIERS:
            raise ValueError(
                f"Unsupported presentation tier for {signal_id}: {tier}"
            )
        matches = reviewed["signal_id"].eq(signal_id)
        if int(matches.sum()) != 1:
            raise ValueError(
                f"Disposition signal must match one evidence row: {signal_id}"
            )
        automated_status = reviewed.loc[
            matches,
            "automated_screen_status",
        ].iloc[0]
        if automated_status not in {
            "candidate_for_human_review",
            "candidate_internal_performance_only",
        }:
            raise ValueError(
                f"Disposition targets a non-candidate signal: {signal_id}"
            )
        reviewed.loc[matches, "human_disposition"] = disposition
        reviewed.loc[matches, "presentation_tier"] = tier
        reviewed.loc[matches, "disposition_rationale"] = str(
            decision["rationale"]
        )
        reviewed.loc[matches, "decision_reviewer"] = str(
            metadata["reviewer"]
        )
        reviewed.loc[matches, "decision_date"] = str(
            metadata["review_date"]
        )

    priority = reviewed["automated_screen_status"].isin(
        [
            "candidate_for_human_review",
            "candidate_internal_performance_only",
        ]
    )
    candidate_dispositions_complete = bool(
        reviewed.loc[priority, "human_disposition"].ne("pending").all()
    )
    metadata = dict(metadata)
    metadata["source_path"] = str(dispositions_path)
    metadata["candidate_dispositions_complete"] = (
        candidate_dispositions_complete
    )
    metadata["recorded_decisions"] = len(decisions)
    metadata["draft_asset_decisions"] = asset_decisions
    metadata["final_approval"] = final_approval
    approved_signal_ids = [
        str(value)
        for value in final_approval.get(
            "approved_public_signal_ids",
            [],
        )
    ]
    for signal_id in approved_signal_ids:
        matches = reviewed["signal_id"].eq(signal_id)
        if int(matches.sum()) != 1:
            raise ValueError(
                f"Approved public signal must match one row: {signal_id}"
            )
        if reviewed.loc[matches, "human_disposition"].iloc[0] != "retain":
            raise ValueError(
                f"Only retained signals may be publicly approved: {signal_id}"
            )
        reviewed.loc[matches, "public_claim_approved"] = True
    metadata["approved_public_signal_ids"] = approved_signal_ids
    return reviewed, metadata


def build_evidence(
    config: dict,
    corrected_out: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    review_config = config.get("result_review", {})
    criteria = {
        "fdr_threshold": float(review_config.get("fdr_threshold", 0.05)),
        "min_biological_replicates": int(
            review_config.get("min_biological_replicates", 20)
        ),
        "require_effect_ci_excludes_zero": bool(
            review_config.get("require_effect_ci_excludes_zero", True)
        ),
        "require_all_sensitivity_scenarios": bool(
            review_config.get(
                "require_all_sensitivity_scenarios",
                True,
            )
        ),
        "require_direction_consistency": bool(
            review_config.get("require_direction_consistency", True)
        ),
        "require_fdr_support_all_scenarios": bool(
            review_config.get(
                "require_fdr_support_all_scenarios",
                True,
            )
        ),
    }
    tables = corrected_out / "tables"
    sensitivity_manifest = (
        corrected_out / "sensitivity_age" / "sensitivity_manifest.csv"
    )
    manifest, sensitivity = _load_sensitivity(sensitivity_manifest)
    expected_scenarios = len(manifest)
    comparison = _load_comparison(
        corrected_out / "comparison" / "association_comparison.csv"
    )
    support = _composition_support(
        tables / "age_celltype_fraction_by_donor.csv"
    )

    rows: list[dict[str, Any]] = []
    rows.extend(
        _association_rows(
            "composition",
            tables / "age_celltype_trend_stats.csv",
            support,
            sensitivity["composition"],
            comparison,
            criteria,
            expected_scenarios,
        )
    )
    rows.extend(
        _association_rows(
            "signature",
            tables / "signature_age_associations.csv",
            support,
            sensitivity["signature"],
            comparison,
            criteria,
            expected_scenarios,
        )
    )
    rows.append(
        _prediction_row(
            tables / "age_pred_metrics.csv",
            corrected_out
            / "comparison"
            / "age_prediction_comparison.csv",
            corrected_out
            / "validation"
            / "replicate_correction_validation.json",
        )
    )
    evidence = pd.DataFrame(rows)
    disposition_metadata: dict[str, Any] = {}
    dispositions_value = review_config.get("dispositions_path")
    if dispositions_value:
        evidence, disposition_metadata = apply_human_dispositions(
            evidence,
            Path(str(dispositions_value)),
        )
    priority = evidence["automated_screen_status"].isin(
        [
            "candidate_for_human_review",
            "candidate_internal_performance_only",
        ]
    )
    candidate_dispositions_complete = bool(
        evidence.loc[priority, "human_disposition"].ne("pending").all()
    )
    public_changes_authorized = bool(
        disposition_metadata.get("status") == "public_promotion_approved"
        and _as_bool(
            disposition_metadata.get("public_claims_approved", False)
        )
        and _as_bool(
            disposition_metadata.get(
                "asset_replacements_approved",
                False,
            )
        )
    )
    summary = {
        "review_status": disposition_metadata.get(
            "status",
            "pending_human_signoff",
        ),
        "public_changes_authorized": public_changes_authorized,
        "corrected_output_root": str(corrected_out),
        "grouping_contract": ASSOCIATION_GROUP,
        "criteria": criteria,
        "sensitivity_scenarios": expected_scenarios,
        "evidence_rows": int(len(evidence)),
        "status_counts": {
            str(key): int(value)
            for key, value in evidence[
                "automated_screen_status"
            ].value_counts().items()
        },
        "analysis_counts": {
            str(key): int(value)
            for key, value in evidence["analysis"].value_counts().items()
        },
        "all_human_dispositions_pending": bool(
            evidence["human_disposition"].eq("pending").all()
        ),
        "candidate_dispositions_complete": (
            candidate_dispositions_complete
        ),
        "human_disposition_counts": {
            str(key): int(value)
            for key, value in evidence.loc[
                priority,
                "human_disposition",
            ].value_counts().items()
        },
        "review_metadata": disposition_metadata,
        "any_public_claim_approved": bool(
            evidence["public_claim_approved"].map(_as_bool).any()
        ),
    }
    return evidence, summary


def _write_json(data: dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    if frame.empty:
        return "_No rows met the automated review-priority criteria._"
    display = frame.loc[:, columns].copy()
    for column in display.select_dtypes(include=["float"]).columns:
        display[column] = display[column].map(
            lambda value: (
                ""
                if pd.isna(value)
                else f"{value:.4g}"
            )
        )
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = [
        "| "
        + " | ".join(str(value).replace("|", "\\|") for value in row)
        + " |"
        for row in display.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def write_review_report(
    evidence: pd.DataFrame,
    summary: dict[str, Any],
    path: Path,
) -> None:
    candidates = evidence.loc[
        evidence["automated_screen_status"].eq(
            "candidate_for_human_review"
        )
    ].sort_values(["analysis", "primary_fdr", "signal_id"])
    prediction = evidence.loc[evidence["analysis"].eq("prediction")].iloc[0]
    criteria = summary["criteria"]
    status_lines = "\n".join(
        f"- `{status}`: {count}"
        for status, count in summary["status_counts"].items()
    )
    candidate_table = _markdown_table(
        candidates,
        [
            "analysis",
            "cell_type",
            "signature",
            "n_replicates",
            "effect_per_10y",
            "effect_ci_low",
            "effect_ci_high",
            "primary_fdr",
            "sensitivity_fdr_supported_count",
            "sensitivity_expected_scenarios",
            "human_disposition",
            "presentation_tier",
        ],
    )
    review_metadata = summary.get("review_metadata", {})
    dispositions_recorded = summary.get(
        "candidate_dispositions_complete",
        False,
    )
    if summary["review_status"] == "public_promotion_approved":
        review_heading = (
            "Corrected public wording and assets approved."
        )
    elif summary["review_status"].startswith(
        "draft_asset_design_approved"
    ):
        review_heading = (
            "Candidate dispositions and draft-asset design approved; "
            "final public sign-off pending."
        )
    elif dispositions_recorded:
        review_heading = (
            "Candidate dispositions approved; public wording and assets "
            "pending."
        )
    else:
        review_heading = "Pending human sign-off."
    reviewer_lines = ""
    if review_metadata:
        reviewer_lines = (
            f'- Reviewer: `{review_metadata["reviewer"]}`\n'
            f'- Review date: `{review_metadata["review_date"]}`\n'
            f'- Public scope: `{review_metadata["public_scope"]}`\n'
        )
    if summary["public_changes_authorized"]:
        review_context = (
            "This document records the pre-specified evidence screen, human "
            "dispositions, and final approval. Public interpretation is "
            "limited to the exact retained signals, wording, and assets "
            "listed in the promotion record."
        )
    else:
        review_context = (
            "This document combines a pre-specified automated evidence "
            "screen with recorded human dispositions. Those dispositions "
            "determine which results may proceed to drafting; they do not "
            "approve biological claims, captions, README text, or public "
            "figure replacement."
        )
    text = f"""# Corrected Results Review

## Review Status

**{review_heading}**

{review_context}

- Corrected output root: `{summary["corrected_output_root"]}`
- Required grouping field: `{summary["grouping_contract"]}`
- Evidence rows: {summary["evidence_rows"]}
{reviewer_lines}- Candidate dispositions complete: **{"Yes" if dispositions_recorded else "No"}**
- Public changes authorized: **{"Yes" if summary["public_changes_authorized"] else "No"}**

## Screening Criteria

Association rows are candidates for human review only when all configured
conditions pass:

- corrected primary FDR < {criteria["fdr_threshold"]};
- corrected effect interval excludes zero;
- at least {criteria["min_biological_replicates"]} biological replicates;
- testable in all {summary["sensitivity_scenarios"]} sensitivity scenarios;
- effect direction agrees in all scenarios;
- primary FDR remains below the threshold in all scenarios;
- the main and sensitivity tables use `biological_replicate_id`.

These criteria are deliberately conservative. Failure means that a row is not
prioritized or remains sensitivity-limited; it does not prove the null.

## Automated Screen

{status_lines}

### Association Candidates

{candidate_table}

`retain` means eligible for cautious drafting, not approved public wording.
`exploratory` results must not be used as headline findings. The excluded
`CD8a/a` result had only 144 cells in total and a median of one cell per
contributing replicate despite passing the automated association screen.

## Internal Age-Prediction Evidence

The selected `{prediction["prediction_model"]}` model used grouped
cross-validation over {int(prediction["n_replicates"])} biological replicates
and {int(prediction["n_splits"])} outer folds. Its donor-level MAE was
{prediction["mae_years"]:.2f} years (95% interval
{prediction["mae_ci_low"]:.2f}-{prediction["mae_ci_high"]:.2f}), versus a
fold-specific mean-age baseline MAE of
{prediction["baseline_mae_years"]:.2f} years; the internal improvement was
{prediction["delta_mae_vs_baseline"]:.2f} years. R-squared was
{prediction["r_squared"]:.3f} (95% interval
{prediction["r2_ci_low"]:.3f}-{prediction["r2_ci_high"]:.3f}).

This supports only an internally cross-validated predictive signal. It is not
external validation and must not be described as a clinical biomarker or
validated aging clock. Its approved presentation tier is
`{prediction["presentation_tier"]}`.

## Interpretation Boundaries

- Composition effects are cross-sectional fraction changes per 10 years, not
  causal effects of aging.
- Signature scores are configured gene-set proxies, not direct measurements of
  pathway activation or suppression.
- FDR, confidence intervals, support, and sensitivity are screening evidence;
  biological plausibility and confounding still require human assessment.
- Provisional-versus-corrected comparison status is context only and is not a
  selection criterion.
- Age, batch, sex, annotation uncertainty, sparse populations, and unequal
  cells per replicate remain possible limitations.

## Human Review Checklist

- [x] Reconcile candidate numbers with the source CSVs.
- [x] Review replicate, age-range, and cell-count support per candidate.
- [x] Assign each prioritized row `retain`, `exploratory`, `exclude`, or
      `defer`.
- [{"x" if summary["public_changes_authorized"] else " "}] Approve exact claim wording and captions without causal or clinical language.
- [{"x" if summary["public_changes_authorized"] else " "}] Approve each proposed asset replacement.
- [x] Record reviewer, date, and rationale in the decision card.

## Promotion State

The machine-readable promotion manifest records the exact authorization state.
The approved wording and captions are preserved in
`docs/reviews/corrected_public_promotion_record.md`; generated analytical
outputs remain separate from curated public assets.
"""
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def build_promotion_manifest(
    evidence: pd.DataFrame,
    corrected_out: Path,
    summary: dict[str, Any],
) -> dict[str, Any]:
    review_metadata = summary.get("review_metadata", {})
    candidates = evidence.loc[
        evidence["automated_screen_status"].isin(
            [
                "candidate_for_human_review",
                "candidate_internal_performance_only",
            ]
        )
    ]
    asset_decisions = review_metadata.get("draft_asset_decisions", {})
    review_figure_dir = corrected_out / "review" / "figures"
    assets = [
        {
            "role": "composition",
            "source": str(
                review_figure_dir / "composition_core_trends_draft.png"
            ),
            "proposed_destination": str(
                Path("docs/assets")
                / "gse164378_corrected_composition_core_trends.png"
            ),
            "draft_decision": asset_decisions.get(
                "composition_figure",
                "",
            ),
            "status": "draft_generation_approved",
        },
        {
            "role": "signature",
            "source": "",
            "proposed_destination": "",
            "draft_decision": asset_decisions.get(
                "signature_figure",
                "",
            ),
            "status": "not_selected_for_curated_assets",
        },
        {
            "role": "prediction",
            "source": str(
                review_figure_dir / "age_prediction_internal_cv_draft.png"
            ),
            "proposed_destination": str(
                Path("docs/assets")
                / "gse164378_corrected_age_prediction_internal_cv.png"
            ),
            "draft_decision": asset_decisions.get(
                "prediction_figure",
                "",
            ),
            "status": "draft_generation_approved",
        },
        {
            "role": "composition_forest",
            "source": str(
                review_figure_dir / "composition_effect_forest_draft.png"
            ),
            "proposed_destination": str(
                Path("docs/assets")
                / "gse164378_corrected_composition_effect_forest.png"
            ),
            "draft_decision": asset_decisions.get(
                "forest_figure",
                "",
            ),
            "status": "draft_generation_approved",
        },
    ]
    retained = evidence.loc[evidence["human_disposition"].eq("retain")]
    exploratory = evidence.loc[
        evidence["human_disposition"].eq("exploratory")
    ]
    excluded = evidence.loc[evidence["human_disposition"].eq("exclude")]
    public_changes_authorized = bool(
        summary.get("public_changes_authorized", False)
    )
    asset_status = (
        "approved_for_public_replacement"
        if public_changes_authorized
        else "draft_generation_approved"
    )
    for asset in assets:
        if asset["role"] != "signature":
            asset["status"] = asset_status
    return {
        "review_status": summary["review_status"],
        "public_changes_authorized": public_changes_authorized,
        "public_scope": review_metadata.get("public_scope", ""),
        "reviewer": review_metadata.get("reviewer", ""),
        "review_date": review_metadata.get("review_date", ""),
        "candidate_signal_ids": candidates["signal_id"].tolist(),
        "retained_for_drafting_signal_ids": retained[
            "signal_id"
        ].tolist(),
        "exploratory_signal_ids": exploratory["signal_id"].tolist(),
        "excluded_signal_ids": excluded["signal_id"].tolist(),
        "assets": assets,
        "required_human_dispositions": [
            "retain",
            "exploratory",
            "exclude",
            "defer",
        ],
        "prohibited_claim_types": [
            "causal aging effect",
            "clinical biomarker",
            "validated aging clock",
            "rejuvenation effect",
            "external validation",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a traceable, non-promotional review of corrected results."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--corrected-out", required=True)
    parser.add_argument("--evidence-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--promotion-out", required=True)
    parser.add_argument("--report-out", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    corrected_out = Path(args.corrected_out)
    evidence, summary = build_evidence(config, corrected_out)

    evidence_path = Path(args.evidence_out)
    ensure_dir(evidence_path.parent)
    evidence.to_csv(evidence_path, index=False)
    _write_json(summary, Path(args.summary_out))
    _write_json(
        build_promotion_manifest(evidence, corrected_out, summary),
        Path(args.promotion_out),
    )
    write_review_report(evidence, summary, Path(args.report_out))


if __name__ == "__main__":
    main()
