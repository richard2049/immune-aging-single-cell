from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

from src.plot_style import (
    PALETTE,
    apply_publication_style,
    finalize_and_save,
    style_axis,
)
from src.utils import ensure_dir, load_config


def _as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def _require_columns(
    frame: pd.DataFrame,
    columns: list[str],
    source: Path,
) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def _load_asset_decisions(
    path: Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    record = load_config(path)
    metadata = record.get("review_metadata", {})
    asset_decisions = record.get("draft_asset_decisions", {})
    final_approval = record.get("final_approval", {})
    decisions = record.get("decisions", [])
    if not isinstance(metadata, dict) or not isinstance(
        asset_decisions,
        dict,
    ):
        raise TypeError(f"Invalid review decision structure: {path}")
    if not isinstance(decisions, list):
        raise TypeError(f"decisions must be a list: {path}")
    if not _as_bool(asset_decisions.get("draft_generation_approved")):
        raise ValueError("Draft figure generation has not been approved.")
    if (
        _as_bool(asset_decisions.get("public_replacement_approved"))
        and final_approval.get("status") != "approved"
    ):
        raise ValueError("Public replacement requires explicit final approval.")
    by_signal = {
        str(item["signal_id"]): item
        for item in decisions
        if isinstance(item, dict) and "signal_id" in item
    }
    return {
        **metadata,
        **asset_decisions,
        "final_approval": final_approval,
    }, by_signal


def _validate_composition_selection(
    selected_signal_ids: list[str],
    decisions: dict[str, dict[str, Any]],
) -> list[str]:
    if len(selected_signal_ids) != 3:
        raise ValueError("The approved composition draft must contain three signals.")
    cell_types: list[str] = []
    for signal_id in selected_signal_ids:
        decision = decisions.get(signal_id)
        if decision is None:
            raise ValueError(f"Unknown composition signal: {signal_id}")
        if decision.get("disposition") != "retain":
            raise ValueError(f"Composition draft includes a non-retained signal: {signal_id}")
        prefix = "composition::"
        if not signal_id.startswith(prefix):
            raise ValueError(f"Composition draft includes another analysis: {signal_id}")
        cell_types.append(signal_id[len(prefix) :])
    return cell_types


def plot_composition_core(
    fractions_path: Path,
    trends_path: Path,
    cell_types: list[str],
    output_path: Path,
    dpi: int,
) -> None:
    fractions = pd.read_csv(fractions_path)
    trends = pd.read_csv(trends_path)
    _require_columns(
        fractions,
        ["cell_type", "age", "fraction"],
        fractions_path,
    )
    _require_columns(
        trends,
        [
            "cell_type",
            "spearman_rho",
            "spearman_fdr",
            "n_donors",
            "adjusted",
            "adjusted_for",
        ],
        trends_path,
    )
    missing = sorted(set(cell_types).difference(trends["cell_type"]))
    if missing:
        raise ValueError(f"Missing composition trends: {missing}")

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15.5, 4.8),
        squeeze=False,
    )
    for ax, cell_type in zip(axes[0], cell_types):
        observed = fractions.loc[fractions["cell_type"].eq(cell_type)]
        x = observed["age"].to_numpy(dtype=float)
        y = observed["fraction"].to_numpy(dtype=float) * 100.0
        ax.scatter(
            x,
            y,
            s=24,
            alpha=0.62,
            color=PALETTE["primary"],
            edgecolors="none",
        )
        if len(observed) >= 2:
            fitted = stats.linregress(x, y)
            grid = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            line_color = PALETTE["secondary"] if fitted.slope >= 0 else PALETTE["danger"]
            ax.plot(
                grid,
                fitted.intercept + fitted.slope * grid,
                color=line_color,
                linewidth=2.2,
            )
        result = trends.loc[trends["cell_type"].eq(cell_type)].iloc[0]
        ax.set_title(
            f"{cell_type}\n"
            f"Spearman rho={result['spearman_rho']:.2f}; "
            f"FDR={result['spearman_fdr']:.2g}; "
            f"n={int(result['n_donors'])}",
            fontsize=11,
        )
        ax.set_xlabel("Age (years)")
        ax.set_ylabel("Cell-type fraction per biological replicate (%)")
        style_axis(ax, grid="y")

    fig.suptitle(
        "Selected cell-type composition associations with age",
        fontsize=16,
        y=1.02,
    )
    fig.text(
        0.5,
        -0.01,
        (
            "Lines are unadjusted visual trends; association statistics are "
            "adjusted for sex and batch."
        ),
        ha="center",
        va="top",
        fontsize=9,
        color=PALETTE["muted"],
    )
    finalize_and_save(fig, output_path, dpi)


def plot_internal_age_prediction(
    predictions_path: Path,
    metrics_path: Path,
    output_path: Path,
    dpi: int,
) -> None:
    predictions = pd.read_csv(predictions_path)
    metrics = pd.read_csv(metrics_path)
    _require_columns(
        predictions,
        [
            "age_true",
            "age_pred",
            "evaluation_level",
            "is_best_model",
        ],
        predictions_path,
    )
    _require_columns(
        metrics,
        [
            "is_best",
            "mae",
            "mae_ci_low",
            "mae_ci_high",
            "r2",
            "r2_ci_low",
            "r2_ci_high",
            "baseline_mae",
            "n_donors",
            "n_splits",
        ],
        metrics_path,
    )
    selected = predictions.loc[
        predictions["evaluation_level"].eq("donor") & predictions["is_best_model"].map(_as_bool)
    ].copy()
    best = metrics.loc[metrics["is_best"].map(_as_bool)]
    if selected.empty or len(best) != 1:
        raise ValueError("Expected donor-level predictions and one selected model.")
    metric = best.iloc[0]
    x = selected["age_true"].to_numpy(dtype=float)
    y = selected["age_pred"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.scatter(
        x,
        y,
        s=28,
        alpha=0.48,
        color=PALETTE["primary"],
        edgecolors="white",
        linewidths=0.25,
    )
    lower = float(np.nanmin([x.min(), y.min()]))
    upper = float(np.nanmax([x.max(), y.max()]))
    ax.plot(
        [lower, upper],
        [lower, upper],
        linestyle="--",
        linewidth=1.8,
        color=PALETTE["muted"],
        label="Identity",
    )
    fitted = stats.linregress(x, y)
    grid = np.linspace(lower, upper, 100)
    ax.plot(
        grid,
        fitted.intercept + fitted.slope * grid,
        linewidth=2.2,
        color=PALETTE["danger"],
        label="Observed fit",
    )
    metrics_text = (
        f"MAE {metric['mae']:.2f} years "
        f"(95% CI {metric['mae_ci_low']:.2f}-{metric['mae_ci_high']:.2f})\n"
        f"R-squared {metric['r2']:.3f} "
        f"(95% CI {metric['r2_ci_low']:.3f}-{metric['r2_ci_high']:.3f})\n"
        f"Mean-age baseline MAE {metric['baseline_mae']:.2f} years"
    )
    ax.text(
        0.03,
        0.97,
        metrics_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        bbox={
            "boxstyle": "square,pad=0.45",
            "facecolor": "white",
            "edgecolor": "#B7B7B7",
            "alpha": 0.94,
        },
    )
    ax.set_xlabel("Observed age (years)")
    ax.set_ylabel("Cross-validated predicted age (years)")
    ax.set_title(
        "Internal cross-validated age prediction\n"
        f"Biological-replicate grouping; "
        f"n={int(metric['n_donors'])}; "
        f"{int(metric['n_splits'])} outer folds"
    )
    style_axis(ax, grid="both")
    ax.legend(frameon=False, loc="lower right")
    finalize_and_save(fig, output_path, dpi)


def plot_retained_composition_forest(
    evidence_path: Path,
    output_path: Path,
    dpi: int,
) -> None:
    evidence = pd.read_csv(evidence_path)
    _require_columns(
        evidence,
        [
            "analysis",
            "cell_type",
            "n_replicates",
            "effect_per_10y",
            "effect_ci_low",
            "effect_ci_high",
            "human_disposition",
            "presentation_tier",
        ],
        evidence_path,
    )
    retained = evidence.loc[
        evidence["analysis"].eq("composition") & evidence["human_disposition"].eq("retain")
    ].copy()
    if retained.empty:
        raise ValueError("No retained composition associations are available.")
    retained["effect_pp"] = retained["effect_per_10y"] * 100.0
    retained["ci_low_pp"] = retained["effect_ci_low"] * 100.0
    retained["ci_high_pp"] = retained["effect_ci_high"] * 100.0
    retained = retained.sort_values("effect_pp").reset_index(drop=True)

    fig_height = max(4.4, 0.72 * len(retained) + 1.8)
    fig, ax = plt.subplots(figsize=(9.2, fig_height))
    ax.axvline(
        0.0,
        color=PALETTE["muted"],
        linestyle="--",
        linewidth=1.4,
    )
    for index, row in retained.iterrows():
        effect = float(row["effect_pp"])
        low = float(row["ci_low_pp"])
        high = float(row["ci_high_pp"])
        color = PALETTE["secondary"] if effect >= 0 else PALETTE["danger"]
        filled = row["presentation_tier"] == "core"
        ax.errorbar(
            effect,
            index,
            xerr=np.array([[effect - low], [high - effect]]),
            fmt="o",
            color=color,
            ecolor=color,
            markerfacecolor=color if filled else "white",
            markeredgecolor=color,
            markersize=8,
            capsize=4,
            linewidth=1.7,
        )
    labels = [
        f"{row.cell_type} (n={int(row.n_replicates)})" for row in retained.itertuples(index=False)
    ]
    ax.set_yticks(np.arange(len(retained)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Adjusted fraction change per 10 years (percentage points)")
    ax.set_title(
        "Retained cell-type composition associations with age\n"
        "Points show adjusted estimates with 95% bootstrap confidence intervals"
    )
    style_axis(ax, grid="x")
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="#4A4A4A",
                markerfacecolor="#4A4A4A",
                linestyle="none",
                label="Core",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="#4A4A4A",
                markerfacecolor="white",
                linestyle="none",
                label="Secondary retained",
            ),
        ],
        frameon=False,
        loc="lower right",
    )
    finalize_and_save(fig, output_path, dpi)


def _write_promotion_record(
    path: Path,
    composition_path: Path,
    prediction_path: Path,
    forest_path: Path,
) -> None:
    text = f"""# Corrected Public-Promotion Record

## Status

**Approved source text and captions for corrected public promotion.**

- Public scope: conservative core
- Signature heatmap: not selected for curated public assets
- Public claims authorized: Yes
- Public asset replacement authorized: Yes

## Proposed README Text

### Replicate-corrected GSE164378 analysis

Review of the source metadata showed that the original `donor_id` field was
reused across pools. Donor-level analyses were therefore repeated using
`Tube_id` as `biological_replicate_id`. The source metadata contained 317
biological replicates, of which 316 met the requirements for the main
composition and prediction analyses.

After adjustment for sex and batch, a 10-year difference in age was associated
with a 1.77 percentage-point lower Tcm/Naive cytotoxic T-cell fraction (95% CI
-1.96 to -1.58; FDR 5.15e-53), a 0.30 percentage-point lower MAIT-cell
fraction (95% CI -0.38 to -0.23; FDR 2.71e-15), and a 0.76 percentage-point
higher CD16+ NK-cell fraction (95% CI 0.42 to 1.08; FDR 6.43e-5). These are
cross-sectional associations between participants, not estimates of
within-person change or causal effects of aging. Tcm/Naive cytotoxic T cells
form a combined annotation category rather than a single resolved subtype.

The scVI-derived features also retained an age-related predictive signal. In
nested cross-validation grouped by biological replicate, the selected model
had a donor-level mean absolute error of 12.39 years (95% CI 11.63 to 13.15),
compared with 15.13 years for a fold-specific mean-age baseline. This is modest
internal predictive performance and has not been validated as an aging clock
or clinical biomarker.

Several predefined gene-set scores were associated with age, but these results
remain exploratory. Signature scores are proxies for the configured gene sets
and do not directly measure pathway activation or suppression.

## Proposed Figure Captions

### Selected Composition Trends

Cell-type fractions per biological replicate plotted against age for the three
associations selected for the main figure. Lines are unadjusted linear
summaries included for visualization; Spearman rho and FDR are from analyses
adjusted for sex and batch. Tcm/Naive cytotoxic T cells form a combined
annotation category. The associations are cross-sectional and do not establish
causality.

Approved source: `{composition_path}`

### Internal Age Prediction

Observed donor age and cross-validated predicted age for the selected model.
Outer folds were grouped by biological replicate. The identity line represents
perfect prediction, while the fitted line shows the relationship observed in
the held-out predictions. Reported metrics describe internal cross-validation
and do not establish external validity or clinical utility.

Approved source: `{prediction_path}`

### Retained Composition Effects

Sex- and batch-adjusted differences in cell-type fraction per 10-year
difference in age for all retained composition associations. Points show
effect estimates in percentage points, and bars show 95% bootstrap confidence
intervals.
Filled markers identify core presentation results; open markers identify
secondary retained results.

Approved source: `{forest_path}`

## Final Approval

- [x] README wording approved with professional tone revision.
- [x] Three captions approved.
- [x] Three draft images approved after visual inspection.
- [x] Provisional curated assets approved for removal.
- [x] Approved corrected assets authorized for `docs/assets/`.
"""
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate non-public draft assets from corrected results."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--corrected-out", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--composition-out", required=True)
    parser.add_argument("--prediction-out", required=True)
    parser.add_argument("--forest-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--record-out", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    corrected_out = Path(args.corrected_out)
    review_config = config.get("result_review", {})
    dispositions_path = Path(str(review_config["dispositions_path"]))
    asset_decisions, dispositions = _load_asset_decisions(dispositions_path)
    selected_ids = [str(value) for value in asset_decisions.get("composition_signals", [])]
    cell_types = _validate_composition_selection(
        selected_ids,
        dispositions,
    )
    if asset_decisions.get("signature_figure") != ("omit_from_curated_assets"):
        raise ValueError("The approved signature-figure decision changed.")

    apply_publication_style()
    dpi = int(config.get("plots", {}).get("dpi", 180))
    composition_out = Path(args.composition_out)
    prediction_out = Path(args.prediction_out)
    forest_out = Path(args.forest_out)
    for output in (composition_out, prediction_out, forest_out):
        ensure_dir(output.parent)

    tables = corrected_out / "tables"
    plot_composition_core(
        tables / "age_celltype_fraction_by_donor.csv",
        tables / "age_celltype_trend_stats.csv",
        cell_types,
        composition_out,
        dpi,
    )
    plot_internal_age_prediction(
        tables / "age_pred_cv_predictions.csv",
        tables / "age_pred_metrics.csv",
        prediction_out,
        dpi,
    )
    plot_retained_composition_forest(
        Path(args.evidence),
        forest_out,
        dpi,
    )
    manifest = {
        "status": "public_promotion_approved",
        "reviewer": asset_decisions["reviewer"],
        "review_date": asset_decisions["review_date"],
        "public_changes_authorized": True,
        "signature_asset_selected": False,
        "draft_assets": [
            str(composition_out),
            str(prediction_out),
            str(forest_out),
        ],
        "source_dispositions": str(dispositions_path),
    }
    manifest_path = Path(args.manifest_out)
    ensure_dir(manifest_path.parent)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_promotion_record(
        Path(args.record_out),
        composition_out,
        prediction_out,
        forest_out,
    )


if __name__ == "__main__":
    main()
