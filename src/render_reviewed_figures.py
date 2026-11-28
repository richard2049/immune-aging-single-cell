"""Render the three reviewed Blood Age Atlas figures from public aggregates."""

from __future__ import annotations

import argparse
import csv
from math import isfinite
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from src.plot_style import PALETTE, apply_publication_style, finalize_and_save, style_axis

ASSETS = Path(__file__).resolve().parents[1] / "docs" / "assets"
DPI = 220
PREDICTION_ROWS = {
    "Nested CV selection": ("visit_mae", "participant bootstrap"),
    "Ridge (α = 0.1)": ("visit_mae", "participant bootstrap"),
    "Age-only baseline": ("visit_mae", "none"),
    "Ridge minus nested CV selection": (
        "paired_fold_mae_difference",
        "paired outer-fold bootstrap",
    ),
}
CD8_ROWS = (
    "Primary (repeated-visits GEE)",
    "Sensitivity (one visit per person)",
)


def read_rows(path, key, expected_count, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter=delimiter))
    keys = [row[key] for row in rows]
    if len(rows) != expected_count or len(set(keys)) != len(rows) or not all(keys):
        raise ValueError(f"{path.name} must contain {expected_count} uniquely labeled rows.")
    return rows


def interval(row, estimate, low, high):
    value, lower, upper = (float(row[field]) for field in (estimate, low, high))
    if not all(map(isfinite, (value, lower, upper))) or not lower <= value <= upper:
        raise ValueError(f"Invalid estimate or interval in {row}.")
    return value, lower, upper


def composition_inputs(path):
    rows = read_rows(path, "cell_type", 10)
    for row in rows:
        if row["association_effect_scale"] != "log_odds_per_10_years":
            raise ValueError("Composition effects must be log-odds per decade.")
    series = []
    for prefix in ("association", "sensitivity"):
        series.append(
            [
                interval(
                    row,
                    f"{prefix}_effect_per_10y",
                    f"{prefix}_effect_ci_low",
                    f"{prefix}_effect_ci_high",
                )
                for row in rows
            ]
        )
    return [row["cell_type"] for row in rows], series


def prediction_inputs(path):
    rows = {row["label"]: row for row in read_rows(path, "label", 4)}
    if rows.keys() != PREDICTION_ROWS.keys():
        raise ValueError("Expected the three visit MAEs and the paired fold comparison.")
    support = {tuple(int(r[k]) for k in ("people", "visits", "outer_folds")) for r in rows.values()}
    if len(support) != 1:
        raise ValueError("Prediction rows disagree on participant, visit, or fold support.")
    people, visits, folds = next(iter(support))
    if not 2 <= folds <= people <= visits:
        raise ValueError("Invalid prediction participant, visit, or fold support.")
    plotted = []
    for label, (summary, method) in PREDICTION_ROWS.items():
        row = rows[label]
        if (row["summary"], row["interval_method"]) != (summary, method):
            raise ValueError(f"Invalid prediction metric or interval method for {label}.")
        if method == "none":
            value = float(row["estimate_years"])
            if not isfinite(value) or value < 0 or row["ci95_low_years"] or row["ci95_high_years"]:
                raise ValueError("The training-age baseline must be a nonnegative point estimate.")
            estimate = value, None, None
        else:
            estimate = interval(row, "estimate_years", "ci95_low_years", "ci95_high_years")
        if summary == "visit_mae":
            if any(value is not None and value < 0 for value in estimate):
                raise ValueError(f"Negative MAE or interval for {label}.")
            plotted.append((label if method != "none" else "Training-age baseline", estimate))
    return plotted, f"{people} people, {visits} visits, {folds} person-grouped outer folds"


def cd8_inputs(path):
    rows = {row["analysis"]: row for row in read_rows(path, "analysis", 2, delimiter="\t")}
    if rows.keys() != set(CD8_ROWS):
        raise ValueError("Expected the primary and one-visit CD8 models.")
    plotted = []
    for label in CD8_ROWS:
        row = rows[label]
        people, visits = int(row["people"]), int(row["visits"])
        if not 0 < people <= visits or (label == CD8_ROWS[1] and visits != people):
            raise ValueError(f"Invalid participant and visit support for {label}.")
        plotted.append(
            (
                f"{label}\n{people} people, {visits} visits",
                interval(row, "estimate", "ci95_low", "ci95_high"),
            )
        )
    return plotted


def plot_interval(ax, y, estimate, marker, color, label=None):
    value, low, high = estimate
    if low is None:
        ax.plot(value, y, marker=marker, color=color)
    else:
        ax.errorbar(
            value,
            y,
            xerr=[[value - low], [high - value]],
            fmt=marker,
            color=color,
            capsize=3,
            linewidth=1.6,
            label=label,
        )


def composition_figure(labels, series):
    fig, ax = plt.subplots(figsize=(8, 5.4))
    for values, offset, marker, color, meaning in zip(
        series,
        (-0.15, 0.15),
        ("o", "s"),
        (PALETTE["primary"], PALETTE["accent"]),
        ("Primary: repeated visits", "Sensitivity: one visit per person"),
        strict=True,
    ):
        for index, estimate in enumerate(values):
            plot_interval(
                ax, index + offset, estimate, marker, color, meaning if index == 0 else None
            )
    ax.set_yticks(range(len(labels)), labels)
    ax.invert_yaxis()
    ax.axvline(0, linestyle="--", color=PALETTE["muted"])
    ax.set_xlabel("Adjusted log-odds per decade of age (95% CI)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1), ncol=2, frameon=False)
    fig.suptitle("Immune-cell composition")
    style_axis(ax, grid="x")
    return fig


def prediction_figure(rows, cohort):
    fig, ax = plt.subplots(figsize=(8, 3.4))
    for index, (label, estimate) in enumerate(rows):
        value, low, high = estimate
        plot_interval(
            ax,
            index,
            estimate,
            "D" if low is None else "o",
            PALETTE["primary"] if index == 0 else PALETTE["muted"],
        )
        display = f"{value:.2f}" if low is None else f"{value:.2f} [{low:.2f}, {high:.2f}]"
        ax.annotate(
            display,
            (value if high is None else high, index),
            xytext=(6, 0),
            textcoords="offset points",
            va="center",
        )
    ax.margins(x=0.24)
    ax.set_yticks(range(len(rows)), [label for label, _ in rows])
    ax.set_ylim(len(rows) - 0.55, -0.45)
    ax.set_xlabel("Visit-level mean absolute error (years)")
    ax.set_title(cohort, fontsize=10, fontweight="normal")
    fig.suptitle("Internal chronological age prediction", weight="semibold")
    style_axis(ax, grid="x")
    return fig


def cd8_figure(rows):
    fig, ax = plt.subplots(figsize=(8, 3.2))
    for index, (_, estimate) in enumerate(rows):
        plot_interval(
            ax, index, estimate, ("o", "s")[index], (PALETTE["primary"], PALETTE["teal"])[index]
        )
    ax.set_yticks(range(len(rows)), [label for label, _ in rows])
    ax.set_ylim(len(rows) - 0.55, -0.45)
    ax.axvline(0, linestyle="--", color=PALETTE["muted"])
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.set_xlabel("Signature-score change per decade of age (95% CI)")
    ax.set_title("Adjusted for sex and batch", fontsize=10, fontweight="normal")
    fig.suptitle("CD8 T-cell SASP-related RNA score", weight="semibold")
    style_axis(ax, grid="x")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True, help="A new output directory.")
    args = parser.parse_args()
    destination = args.out_dir.resolve()
    if destination == ASSETS or ASSETS in destination.parents:
        parser.error("Choose an output directory outside the approved asset directory.")
    try:
        composition = composition_inputs(ASSETS / "blood_age_atlas_composition_summary.csv")
        prediction = prediction_inputs(ASSETS / "blood_age_atlas_age_prediction_summary.csv")
        cd8 = cd8_inputs(ASSETS / "blood_age_atlas_cd8_sasp_summary.tsv")
        destination.mkdir(parents=True, exist_ok=False)
    except (KeyError, TypeError, ValueError, OSError) as exc:
        parser.error(str(exc))
    apply_publication_style(DPI)
    for name, plot, inputs in (
        ("composition_effects", composition_figure, composition),
        ("age_prediction", prediction_figure, prediction),
        ("cd8_sasp", cd8_figure, (cd8,)),
    ):
        finalize_and_save(plot(*inputs), destination / f"blood_age_atlas_{name}.png", DPI)


if __name__ == "__main__":
    main()
