from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import scanpy as sc

from .plot_style import PALETTE, apply_publication_style, finalize_and_save, save_placeholder, style_axis
from .utils import ensure_dir, load_config


def _save_placeholder(path: Path, title: str, message: str, dpi: int) -> None:
    save_placeholder(path=path, title=title, message=message, dpi=dpi)


def _select_label_column(adata, requested: str = "auto") -> str | None:
    if requested and requested.lower() != "auto":
        return requested if requested in adata.obs.columns else None

    candidates = ["cell_type", "majority_voting", "predicted_labels", "Cluster_names", "leiden"]
    for c in candidates:
        if c in adata.obs.columns:
            return c
    return None


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
        size=12,
    )
    fig = plt.gcf()
    if fig.axes:
        fig.axes[0].set_title(f"UMAP colored by {color} (n={adata.n_obs:,} cells)")
    finalize_and_save(fig, path, dpi)


def _save_qc_distributions(adata, path: Path, dpi: int) -> None:
    qc_candidates = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "percent.ribo",
    ]
    qc_cols = [c for c in qc_candidates if c in adata.obs.columns]
    qc_cols = list(dict.fromkeys(qc_cols))[:4]

    if not qc_cols:
        _save_placeholder(path, "QC Distributions", "No QC columns found in adata.obs.", dpi)
        return

    n_cols = 2
    n_rows = int(np.ceil(len(qc_cols) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 4 * n_rows))
    axes = np.array(axes).reshape(-1)

    for ax, col in zip(axes, qc_cols):
        vals = adata.obs[col].dropna().to_numpy()
        ax.hist(vals, bins=40, color=PALETTE["primary"], edgecolor="white", alpha=0.85)
        median = float(np.median(vals)) if vals.size > 0 else np.nan
        ax.axvline(
            median,
            color=PALETTE["danger"],
            linestyle="--",
            linewidth=1.5,
            label=f"median={median:.2f}",
        )
        ax.set_title(f"{col} (n={vals.size:,})")
        ax.set_xlabel(col)
        ax.set_ylabel("Cells")
        if col in {"total_counts", "nCount_RNA"}:
            ax.set_xscale("log")
            ax.set_xlabel(f"{col} (log scale)")
        style_axis(ax, grid="y")
        ax.legend(frameon=False, fontsize=8)

    for ax in axes[len(qc_cols) :]:
        ax.axis("off")

    finalize_and_save(fig, path, dpi)


def _save_celltype_composition(adata, path: Path, top_n: int, dpi: int, label_col: str | None) -> None:
    if label_col is None:
        _save_placeholder(path, "Cell-type Composition", "Missing label column in adata.obs.", dpi)
        return

    counts = adata.obs[label_col].astype(str).value_counts()
    if counts.empty:
        _save_placeholder(path, "Cell-type Composition", f"{label_col} column is empty.", dpi)
        return

    plot_df = counts.head(top_n).sort_values(ascending=True)
    frac = plot_df / float(adata.n_obs)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(plot_df.index, frac.values, color=PALETTE["teal"])
    ax.set_xlabel("Fraction of cells")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax.set_ylabel("Label")
    title = f"Top {len(plot_df)} labels"
    if label_col != "cell_type":
        title += f" (source: {label_col})"
    ax.set_title(f"{title}\nTotal cells: {adata.n_obs:,}")
    style_axis(ax, grid="x")
    for i, (lbl, val) in enumerate(plot_df.items()):
        ax.text(frac.loc[lbl], i, f"  {int(val):,}", va="center", fontsize=9)
    finalize_and_save(fig, path, dpi)


def _write_tables(adata, table_dir: Path, label_col: str | None) -> None:
    if label_col is not None:
        counts = adata.obs[label_col].astype(str).value_counts().rename_axis("cell_type")
    else:
        counts = pd.Series({"unlabeled": int(adata.n_obs)}, name="n_cells")
        counts.index.name = "cell_type"

    counts_df = counts.rename("n_cells").reset_index()
    counts_df["fraction"] = counts_df["n_cells"] / float(max(adata.n_obs, 1))
    counts_df["fraction_percent"] = counts_df["fraction"] * 100.0
    counts_df["rank"] = np.arange(1, len(counts_df) + 1)
    counts_df["label_source"] = label_col if label_col is not None else "none"
    counts_df.to_csv(table_dir / "cell_type_counts.csv", index=False)

    frac_df = counts_df[["cell_type", "n_cells", "fraction", "fraction_percent", "rank", "label_source"]].copy()
    frac_df.to_csv(table_dir / "cell_type_fractions.csv", index=False)

    if "leiden" in adata.obs and label_col is not None and label_col != "leiden":
        ct = pd.crosstab(
            adata.obs["leiden"].astype(str),
            adata.obs[label_col].astype(str),
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
    label_requested = str(report_cfg.get("label_obs", "auto"))
    use_backed = bool(report_cfg.get("use_backed", True))

    fig_dir = Path(args.figdir)
    table_dir = Path(args.tabledir)
    ensure_dir(fig_dir)
    ensure_dir(table_dir)
    apply_publication_style(dpi=dpi)

    adata = sc.read_h5ad(args.inp, backed="r" if use_backed else None)
    label_col = _select_label_column(adata, requested=label_requested)
    umap_label = label_col if label_col is not None else "cell_type"

    _save_umap(adata, "leiden", fig_dir / "umap_leiden.png", dpi)
    _save_umap(adata, umap_label, fig_dir / "umap_cell_type.png", dpi)
    _save_qc_distributions(adata, fig_dir / "qc_distributions.png", dpi)
    _save_celltype_composition(adata, fig_dir / "cell_type_composition.png", top_n, dpi, label_col=label_col)

    _write_tables(adata, table_dir, label_col=label_col)
    if hasattr(adata, "file") and getattr(adata, "file", None) is not None:
        adata.file.close()

    Path(args.done).write_text("report_complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()
