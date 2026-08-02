from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_dir


def _require_columns(
    table: pd.DataFrame,
    required: set[str],
    label: str,
) -> None:
    missing = sorted(required.difference(table.columns))
    if missing:
        raise KeyError(f"{label} is missing columns: {missing}")


def _check_grouping_provenance(
    table: pd.DataFrame,
    label: str,
    expected: str,
) -> None:
    _require_columns(table, {"grouping_id_column"}, label)
    values = set(table["grouping_id_column"].dropna().astype(str))
    if values != {expected}:
        raise ValueError(
            f"{label} grouping provenance is {sorted(values)}, "
            f"expected only {expected!r}."
        )


def validate(out_dir: Path) -> dict[str, Any]:
    tables = out_dir / "tables"
    with (tables / "biological_replicate_audit.json").open(
        "r",
        encoding="utf-8",
    ) as handle:
        replicate_audit = json.load(handle)
    if not bool(replicate_audit.get("passed")):
        raise ValueError("Biological-replicate audit did not pass.")
    if int(replicate_audit.get("mapped_cells", 0)) != 1_000_000:
        raise ValueError("Expected exactly 1,000,000 mapped cells.")
    if int(replicate_audit.get("unique_biological_replicates", 0)) != 317:
        raise ValueError("Expected exactly 317 biological replicates.")

    mapping = pd.read_csv(
        tables / "cell_to_biological_replicate.csv",
        usecols=["cell_id", "biological_replicate_id"],
        dtype="string",
    )
    if len(mapping) != 1_000_000:
        raise ValueError(f"Mapping has {len(mapping):,} rows, not 1,000,000.")
    if mapping["cell_id"].duplicated().any():
        raise ValueError("Mapping contains duplicate cell IDs.")
    if mapping["biological_replicate_id"].isna().any():
        raise ValueError("Mapping contains missing biological replicate IDs.")

    composition = pd.read_csv(tables / "age_celltype_trend_stats.csv")
    signatures = pd.read_csv(tables / "signature_age_associations.csv")
    prediction_metrics = pd.read_csv(tables / "age_pred_metrics.csv")
    prediction_summary = pd.read_csv(
        tables / "age_pred_model_comparison_summary.csv"
    )
    for label, table in [
        ("composition trends", composition),
        ("signature associations", signatures),
        ("age-prediction metrics", prediction_metrics),
        ("age-prediction summary", prediction_summary),
    ]:
        if table.empty:
            raise ValueError(f"{label} is empty.")
        _check_grouping_provenance(
            table,
            label,
            "biological_replicate_id",
        )

    donor_fractions = pd.read_csv(
        tables / "age_celltype_fraction_by_donor.csv",
        usecols=["biological_replicate_id"],
    )
    signature_scores = pd.read_csv(
        tables / "signature_scores_by_donor_celltype.csv",
        usecols=["biological_replicate_id"],
    )
    predictions = pd.read_csv(
        tables / "age_pred_cv_predictions.csv",
        usecols=[
            "model",
            "params",
            "evaluation_level",
            "biological_replicate_id",
            "fold",
        ],
    )
    fold_counts = predictions.groupby(
        [
            "model",
            "params",
            "evaluation_level",
            "biological_replicate_id",
        ],
        observed=True,
    )["fold"].nunique()
    leakage_groups = int((fold_counts > 1).sum())
    if leakage_groups:
        raise ValueError(
            f"Detected CV leakage in {leakage_groups} model-replicate groups."
        )

    sensitivity_manifest = pd.read_csv(
        out_dir / "sensitivity_age" / "sensitivity_manifest.csv"
    )
    if sensitivity_manifest.empty:
        raise ValueError("Sensitivity manifest is empty.")
    for _, row in sensitivity_manifest.iterrows():
        scenario_composition = pd.read_csv(
            Path(str(row["composition_trends_table"]))
        )
        scenario_signature = pd.read_csv(
            Path(str(row["signature_assoc_table"]))
        )
        _check_grouping_provenance(
            scenario_composition,
            f"scenario {row['scenario']} composition",
            "biological_replicate_id",
        )
        _check_grouping_provenance(
            scenario_signature,
            f"scenario {row['scenario']} signature",
            "biological_replicate_id",
        )

    expected_figures = [
        "age_celltype_composition_by_bin.png",
        "age_celltype_top_trends.png",
        "signature_age_heatmap.png",
        "signature_age_top_associations.png",
        "age_pred_observed_vs_predicted.png",
        "age_pred_mae_by_celltype.png",
        "supp_age_effect_ci_forest.png",
        "supp_age_sensitivity_stability.png",
    ]
    missing_figures = [
        name
        for name in expected_figures
        if not (out_dir / "figures" / name).is_file()
        or (out_dir / "figures" / name).stat().st_size == 0
    ]
    if missing_figures:
        raise FileNotFoundError(
            f"Missing or empty corrected figures: {missing_figures}"
        )

    comparison_summary = out_dir / "comparison" / "summary.json"
    if not comparison_summary.is_file():
        raise FileNotFoundError(
            f"Missing descriptive comparison: {comparison_summary}"
        )

    return {
        "passed": True,
        "scope": "technical validation; no biological interpretation",
        "checkpoint_strategy": (
            "read-only reuse of the existing annotated H5AD; no remediated "
            "expression checkpoint was written"
        ),
        "mapped_cells": int(len(mapping)),
        "unique_biological_replicates": int(
            mapping["biological_replicate_id"].nunique()
        ),
        "composition_tests": int(len(composition)),
        "signature_tests": int(len(signatures)),
        "age_prediction_rows": int(len(predictions)),
        "age_prediction_replicates": int(
            predictions["biological_replicate_id"].nunique()
        ),
        "cv_leakage_groups": leakage_groups,
        "sensitivity_scenarios": int(len(sensitivity_manifest)),
        "nonempty_figures": len(expected_figures),
        "donor_fraction_replicates": int(
            donor_fractions["biological_replicate_id"].nunique()
        ),
        "signature_score_replicates": int(
            signature_scores["biological_replicate_id"].nunique()
        ),
        "biological_interpretation_status": (
            "provisional_pending_human_review"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    report = validate(Path(args.outdir))
    report_path = Path(args.report)
    ensure_dir(report_path.parent)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        "[validate_replicate_correction] "
        f"passed for {report['mapped_cells']:,} cells, "
        f"{report['unique_biological_replicates']:,} replicates, "
        f"{report['cv_leakage_groups']} CV leakage groups",
        flush=True,
    )


if __name__ == "__main__":
    main()
