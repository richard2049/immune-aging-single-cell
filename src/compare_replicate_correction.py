from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import ensure_dir


def _as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def _compare_associations(
    provisional: Path,
    corrected: Path,
    analysis: str,
    keys: list[str],
    effect_col: str,
    fdr_col: str,
) -> pd.DataFrame:
    old = pd.read_csv(provisional)
    new = pd.read_csv(corrected)
    required = set(keys + [effect_col, fdr_col])
    for label, table in [("provisional", old), ("corrected", new)]:
        missing = sorted(required.difference(table.columns))
        if missing:
            raise KeyError(f"{label} {analysis} table is missing columns: {missing}")

    selected = keys + [effect_col, fdr_col]
    old = old[selected].rename(
        columns={
            effect_col: "provisional_effect",
            fdr_col: "provisional_fdr",
        }
    )
    new = new[selected].rename(
        columns={
            effect_col: "corrected_effect",
            fdr_col: "corrected_fdr",
        }
    )
    merged = old.merge(new, on=keys, how="outer", validate="one_to_one")
    merged.insert(0, "analysis", analysis)
    merged["provisional_supported"] = merged["provisional_fdr"].lt(0.05)
    merged["corrected_supported"] = merged["corrected_fdr"].lt(0.05)
    old_sign = np.sign(merged["provisional_effect"])
    new_sign = np.sign(merged["corrected_effect"])
    both = merged["provisional_effect"].notna() & merged["corrected_effect"].notna()
    merged["comparison_status"] = "stable_direction"
    merged.loc[merged["provisional_effect"].isna(), "comparison_status"] = "newly_testable"
    merged.loc[merged["corrected_effect"].isna(), "comparison_status"] = "no_longer_testable"
    merged.loc[both & old_sign.ne(new_sign), "comparison_status"] = "direction_changed"
    merged.loc[
        both & ~merged["provisional_supported"] & merged["corrected_supported"],
        "comparison_status",
    ] = "newly_supported"
    merged.loc[
        both & merged["provisional_supported"] & ~merged["corrected_supported"],
        "comparison_status",
    ] = "lost_support"
    merged["effect_change"] = merged["corrected_effect"] - merged["provisional_effect"]
    return merged


def _compare_prediction_metrics(
    provisional: Path,
    corrected: Path,
) -> pd.DataFrame:
    old = pd.read_csv(provisional)
    new = pd.read_csv(corrected)
    keys = ["model", "params", "primary_level"]
    metrics = ["mae", "rmse", "r2", "spearman_rho", "n_donors", "n_groups"]
    for label, table in [("provisional", old), ("corrected", new)]:
        missing = sorted(set(keys + metrics).difference(table.columns))
        if missing:
            raise KeyError(f"{label} age-prediction metrics are missing columns: {missing}")
    old = old[keys + metrics].add_prefix("provisional_")
    old = old.rename(columns={f"provisional_{key}": key for key in keys})
    new = new[keys + metrics].add_prefix("corrected_")
    new = new.rename(columns={f"corrected_{key}": key for key in keys})
    merged = old.merge(new, on=keys, how="outer", validate="one_to_one")
    merged.insert(0, "analysis", "age_prediction")
    merged["comparison_status"] = np.where(
        merged["provisional_mae"].notna() & merged["corrected_mae"].notna(),
        "comparable",
        "model_set_changed",
    )
    merged["mae_change"] = merged["corrected_mae"] - merged["provisional_mae"]
    return merged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provisional-out", required=True)
    parser.add_argument("--corrected-out", required=True)
    parser.add_argument("--association-out", required=True)
    parser.add_argument("--prediction-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    old = Path(args.provisional_out)
    new = Path(args.corrected_out)
    composition = _compare_associations(
        old / "tables" / "age_celltype_trend_stats.csv",
        new / "tables" / "age_celltype_trend_stats.csv",
        analysis="composition",
        keys=["cell_type"],
        effect_col="slope_per_10y",
        fdr_col="spearman_fdr",
    )
    signatures = _compare_associations(
        old / "tables" / "signature_age_associations.csv",
        new / "tables" / "signature_age_associations.csv",
        analysis="signature",
        keys=["cell_type", "signature"],
        effect_col="effect_per_10y",
        fdr_col="fdr",
    )
    associations = pd.concat(
        [composition, signatures],
        ignore_index=True,
        sort=False,
    )
    prediction = _compare_prediction_metrics(
        old / "tables" / "age_pred_metrics.csv",
        new / "tables" / "age_pred_metrics.csv",
    )

    association_out = Path(args.association_out)
    prediction_out = Path(args.prediction_out)
    summary_out = Path(args.summary_out)
    for path in [association_out, prediction_out, summary_out]:
        ensure_dir(path.parent)
    associations.to_csv(association_out, index=False)
    prediction.to_csv(prediction_out, index=False)

    status_counts = (
        associations.groupby(["analysis", "comparison_status"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
        .to_dict(orient="records")
    )
    summary = {
        "scope": "descriptive comparison only; no biological interpretation",
        "provisional_out_dir": str(old),
        "corrected_out_dir": str(new),
        "association_status_counts": status_counts,
        "age_prediction_models_compared": int(
            prediction["comparison_status"].eq("comparable").sum()
        ),
        "corrected_biological_interpretation_status": ("provisional_pending_human_review"),
    }
    with summary_out.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
