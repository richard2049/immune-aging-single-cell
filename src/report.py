from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from .utils import ensure_dir, load_config


def _save_placeholder(path: Path, title: str, message: str, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=14, weight="bold")
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_umap(adata, color: str, path: Path, dpi: int) -> None:
    if "X_umap" not in adata.obsm:
        _save_placeholder(path, f"UMAP colored by {color}", "Missing X_umap embedding.", dpi)
        return
    if color not in adata.obs:
        _save_placeholder(path, f"UMAP colored by {color}", f"Missing obs column: {color}", dpi)
        return

    sc.pl.umap(
        adata,
        color=color,
        show=False,
        frameon=False,
        legend_loc="right margin",
    )
    fig = plt.gcf()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_qc_distributions(adata, path: Path, dpi: int) -> None:
    qc_cols = [
        c
        for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
        if c in adata.obs.columns
    ]

    if not qc_cols:
        _save_placeholder(path, "QC Distributions", "No QC columns found in adata.obs.", dpi)
        return

    n_cols = 2
    n_rows = int(np.ceil(len(qc_cols) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 4 * n_rows))
    axes = np.array(axes).reshape(-1)

    for ax, col in zip(axes, qc_cols):
        vals = adata.obs[col].dropna().to_numpy()
        ax.hist(vals, bins=40, color="#4C78A8", edgecolor="white")
        ax.set_title(col)
        ax.set_xlabel(col)
        ax.set_ylabel("Cells")

    for ax in axes[len(qc_cols) :]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_celltype_composition(adata, path: Path, top_n: int, dpi: int) -> None:
    if "cell_type" not in adata.obs:
        _save_placeholder(path, "Cell-type Composition", "Missing obs column: cell_type", dpi)
        return

    counts = adata.obs["cell_type"].astype(str).value_counts()
    if counts.empty:
        _save_placeholder(path, "Cell-type Composition", "cell_type column is empty.", dpi)
        return

    plot_df = counts.head(top_n).sort_values(ascending=True)
    frac = plot_df / float(adata.n_obs)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(plot_df.index, frac.values, color="#72B7B2")
    ax.set_xlabel("Fraction of cells")
    ax.set_ylabel("Cell type")
    ax.set_title(f"Top {len(plot_df)} cell types")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _write_tables(adata, table_dir: Path) -> None:
    if "cell_type" in adata.obs:
        counts = adata.obs["cell_type"].astype(str).value_counts().rename_axis("cell_type")
    else:
        counts = pd.Series({"unlabeled": int(adata.n_obs)}, name="n_cells")
        counts.index.name = "cell_type"

    counts_df = counts.rename("n_cells").reset_index()
    counts_df.to_csv(table_dir / "cell_type_counts.csv", index=False)

    frac_df = counts_df.copy()
    frac_df["fraction"] = frac_df["n_cells"] / float(max(adata.n_obs, 1))
    frac_df.to_csv(table_dir / "cell_type_fractions.csv", index=False)

    if "leiden" in adata.obs and "cell_type" in adata.obs:
        ct = pd.crosstab(
            adata.obs["leiden"].astype(str),
            adata.obs["cell_type"].astype(str),
            rownames=["leiden"],
            colnames=["cell_type"],
        )
    elif "leiden" in adata.obs:
        ct = adata.obs["leiden"].astype(str).value_counts().to_frame(name="n_cells")
        ct.index.name = "leiden"
    else:
        ct = pd.DataFrame({"n_cells": [int(adata.n_obs)]}, index=["all_cells"])
        ct.index.name = "group"

    ct.to_csv(table_dir / "cluster_celltype_crosstab.csv")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabledir", required=True)
    ap.add_argument("--done", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    report_cfg = cfg.get("report", {})

    dpi = int(report_cfg.get("dpi", 150))
    top_n = int(report_cfg.get("max_celltypes_plot", 20))

    fig_dir = Path(args.figdir)
    table_dir = Path(args.tabledir)
    ensure_dir(fig_dir)
    ensure_dir(table_dir)

    adata = sc.read_h5ad(args.inp)

    _save_umap(adata, "leiden", fig_dir / "umap_leiden.png", dpi)
    _save_umap(adata, "cell_type", fig_dir / "umap_cell_type.png", dpi)
    _save_qc_distributions(adata, fig_dir / "qc_distributions.png", dpi)
    _save_celltype_composition(adata, fig_dir / "cell_type_composition.png", top_n, dpi)

    _write_tables(adata, table_dir)

    Path(args.done).write_text("report_complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()
