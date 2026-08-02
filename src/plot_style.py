from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

# Colorblind-safe palette inspired by scientific plotting defaults.
COLOR_CYCLE = [
    "#1F77B4",  # blue
    "#FF7F0E",  # orange
    "#2CA02C",  # green
    "#D62728",  # red
    "#9467BD",  # purple
    "#8C564B",  # brown
    "#E377C2",  # pink
    "#7F7F7F",  # gray
    "#BCBD22",  # olive
    "#17BECF",  # cyan
]

PALETTE = {
    "primary": "#1F77B4",
    "secondary": "#2CA02C",
    "accent": "#FF7F0E",
    "danger": "#D62728",
    "muted": "#6E6E6E",
    "teal": "#17BECF",
    "grid": "#D9D9D9",
}

DIVERGING_CMAP = "RdBu_r"

# Distinct categorical palette for many labels (colorblind-aware, high-contrast order).
DISTINCT_CATEGORICAL_COLORS = [
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
    "#BCBD22",
    "#17BECF",
    "#393B79",
    "#637939",
    "#8C6D31",
    "#843C39",
    "#7B4173",
    "#3182BD",
    "#31A354",
    "#756BB1",
    "#E6550D",
    "#969696",
    "#DD1C77",
    "#6BAED6",
    "#74C476",
    "#FD8D3C",
    "#9E9AC8",
    "#BDBDBD",
]


def apply_publication_style(dpi: int = 180) -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#2B2B2B",
            "axes.linewidth": 0.8,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "axes.titleweight": "semibold",
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "legend.title_fontsize": 9,
            "grid.color": PALETTE["grid"],
            "grid.linestyle": "-",
            "grid.linewidth": 0.6,
            "axes.prop_cycle": cycler("color", COLOR_CYCLE),
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(ax, grid: str | None = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid in {"x", "y", "both"}:
        ax.grid(axis=grid, alpha=0.7)


def finalize_and_save(fig, path: Path, dpi: int) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_placeholder(path: Path, title: str, message: str, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=14, weight="semibold")
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11)
    finalize_and_save(fig, path, dpi)


def categorical_palette(labels: list[str], other_label: str = "Other", other_color: str = "#B0B0B0") -> list[str]:
    labels = [str(x) for x in labels]
    non_other = [x for x in labels if x != other_label]
    colors: dict[str, str] = {}
    for i, lab in enumerate(non_other):
        colors[lab] = DISTINCT_CATEGORICAL_COLORS[i % len(DISTINCT_CATEGORICAL_COLORS)]
    if other_label in labels:
        colors[other_label] = other_color
    return [colors.get(lab, DISTINCT_CATEGORICAL_COLORS[i % len(DISTINCT_CATEGORICAL_COLORS)]) for i, lab in enumerate(labels)]
