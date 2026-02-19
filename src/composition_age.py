from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .plot_style import PALETTE, apply_publication_style, finalize_and_save, save_placeholder, style_axis
from .utils import ensure_dir, load_config


def _save_placeholder(path: Path, title: str, message: str, dpi: int) -> None:
    save_placeholder(path=path, title=title, message=message, dpi=dpi)


def _first_existing(columns: list[str], candidates: list[str]) -> str | None:
    cols = set(columns)
    for key in candidates:
        if key in cols:
            return key
    return None


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(pvals.shape[0], np.nan, dtype=float)
    valid = np.isfinite(pvals)
    if not np.any(valid):
        return out

    p = pvals[valid]
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * float(len(ranked)) / (np.arange(len(ranked)) + 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    restored = np.empty_like(q)
    restored[order] = q
    out[valid] = restored
    return out


def _prepare_obs(inp_h5ad: str, cfg: dict) -> tuple[pd.DataFrame, dict]:
    ccfg = cfg.get("composition_age", {})
    adata = ad.read_h5ad(inp_h5ad, backed="r")
    obs = adata.obs.copy()
    adata.file.close()

    all_cols = list(obs.columns)
    cfg_age_col = ccfg.get("age_col")
    cfg_donor_col = ccfg.get("donor_col")
    cfg_celltype_col = ccfg.get("celltype_col")

    age_col = cfg_age_col if cfg_age_col in obs.columns else _first_existing(all_cols, ["age"])
    donor_col = (
        cfg_donor_col
        if cfg_donor_col in obs.columns
        else _first_existing(all_cols, ["donor_id", "donor_id_y", "donor_id_x", "donor", "sample_id"])
    )
    celltype_col = (
        cfg_celltype_col
        if cfg_celltype_col in obs.columns
        else _first_existing(all_cols, ["cell_type", "majority_voting", "predicted_labels", "leiden"])
    )

    required = {
        "age_col": age_col,
        "donor_col": donor_col,
        "celltype_col": celltype_col,
    }
    if any(v is None for v in required.values()):
        return pd.DataFrame(), required

    out = obs[[age_col, donor_col, celltype_col]].copy()
    out.columns = ["age", "donor_id", "cell_type"]
    out["age"] = pd.to_numeric(out["age"], errors="coerce")
    out["donor_id"] = out["donor_id"].astype(str)
    out["cell_type"] = out["cell_type"].astype(str)
    out = out.dropna(subset=["age", "donor_id", "cell_type"])
    return out, required


def _build_donor_fraction_table(obs: pd.DataFrame, min_cells_per_donor: int) -> pd.DataFrame:
    totals = obs.groupby("donor_id", observed=False).agg(
        age=("age", "median"),
        total_cells=("cell_type", "size"),
    )
    totals = totals[totals["total_cells"] >= min_cells_per_donor].reset_index()
    if totals.empty:
        return pd.DataFrame(columns=["donor_id", "age", "cell_type", "n_cells", "total_cells", "fraction"])

    kept = obs[obs["donor_id"].isin(totals["donor_id"])]
    counts = (
        kept.groupby(["donor_id", "cell_type"], observed=False)
        .size()
        .reset_index(name="n_cells")
    )
    out = counts.merge(totals, on="donor_id", how="left")
    out["fraction"] = out["n_cells"] / out["total_cells"]
    return out[["donor_id", "age", "cell_type", "n_cells", "total_cells", "fraction"]]


def _compute_trends(donor_fraction: pd.DataFrame, min_donors_per_celltype: int) -> pd.DataFrame:
    rows: list[dict] = []
    for ct, g in donor_fraction.groupby("cell_type", observed=False):
        x = g["age"].to_numpy(dtype=float)
        y = g["fraction"].to_numpy(dtype=float)
        if len(g) < min_donors_per_celltype:
            continue
        rho, p = stats.spearmanr(x, y)
        lr = stats.linregress(x, y)
        rows.append(
            {
                "cell_type": ct,
                "n_donors": int(len(g)),
                "mean_fraction": float(np.mean(y)),
                "spearman_rho": float(rho),
                "spearman_pvalue": float(p),
                "slope_per_year": float(lr.slope),
                "slope_per_10y": float(lr.slope * 10.0),
                "slope_pvalue": float(lr.pvalue),
                "r_squared": float(lr.rvalue**2),
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "cell_type",
                "n_donors",
                "mean_fraction",
                "spearman_rho",
                "spearman_pvalue",
                "spearman_fdr",
                "slope_per_year",
                "slope_per_10y",
                "slope_pvalue",
                "slope_fdr",
                "r_squared",
                "direction",
                "fdr_significant",
            ]
        )
    out = pd.DataFrame(rows)
    out["spearman_fdr"] = _bh_fdr(out["spearman_pvalue"].to_numpy())
    out["slope_fdr"] = _bh_fdr(out["slope_pvalue"].to_numpy())
    out["direction"] = np.where(out["slope_per_year"] > 0, "increase_with_age", "decrease_with_age")
    out["fdr_significant"] = out["spearman_fdr"] < 0.05
    return out.sort_values(["spearman_fdr", "spearman_pvalue", "spearman_rho"], ascending=[True, True, False])


def _plot_fraction_by_age_bin(
    donor_fraction: pd.DataFrame,
    path: Path,
    age_bins: list[float],
    top_n_celltypes: int,
    dpi: int,
) -> None:
    if donor_fraction.empty:
        _save_placeholder(path, "Cell-type Composition by Age Bin", "No donor-level fractions available.", dpi)
        return

    df = donor_fraction.copy()
    bins = sorted(set(float(v) for v in age_bins))
    if len(bins) < 2:
        _save_placeholder(path, "Cell-type Composition by Age Bin", "Need at least two age bin edges.", dpi)
        return

    labels = [f"{int(bins[i])}-{int(bins[i + 1])}" for i in range(len(bins) - 1)]
    df["age_bin"] = pd.cut(df["age"], bins=bins, labels=labels, include_lowest=True, right=False)
    df = df.dropna(subset=["age_bin"])
    if df.empty:
        _save_placeholder(path, "Cell-type Composition by Age Bin", "No donors fall inside configured age bins.", dpi)
        return

    top_ct = (
        df.groupby("cell_type", observed=False)["fraction"]
        .mean()
        .sort_values(ascending=False)
        .head(top_n_celltypes)
        .index
    )
    df["cell_type_plot"] = np.where(df["cell_type"].isin(top_ct), df["cell_type"], "Other")

    pivot = (
        df.groupby(["age_bin", "cell_type_plot"], observed=False)["fraction"]
        .mean()
        .unstack(fill_value=0.0)
    )
    pivot = pivot.div(pivot.sum(axis=1).replace(0.0, np.nan), axis=0).fillna(0.0)
    donors_per_bin = df.groupby("age_bin", observed=False)["donor_id"].nunique().reindex(pivot.index).fillna(0).astype(int)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(pivot.index))
    bottom = np.zeros(len(pivot.index), dtype=float)
    for col in pivot.columns:
        y = (pivot[col] * 100.0).to_numpy(dtype=float)
        ax.bar(x, y, bottom=bottom, label=col, width=0.85)
        bottom += y

    ax.set_xticks(x)
    xt = [f"{str(b)}\n(n={int(n)})" for b, n in zip(pivot.index, donors_per_bin.to_numpy())]
    ax.set_xticklabels(xt, rotation=30, ha="right")
    ax.set_ylabel("Mean fraction across donors (%)")
    ax.set_xlabel("Age bin")
    ax.set_ylim(0, 100)
    ax.set_title("Cell-type composition across age bins\nDonor-level mean fractions")
    style_axis(ax, grid="y")
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    finalize_and_save(fig, path, dpi)


def _plot_top_trends(
    donor_fraction: pd.DataFrame,
    trends: pd.DataFrame,
    path: Path,
    top_n_trends: int,
    dpi: int,
) -> None:
    if donor_fraction.empty or trends.empty:
        _save_placeholder(path, "Top Cell-type Age Trends", "No trend statistics available.", dpi)
        return

    top = trends.head(top_n_trends)["cell_type"].tolist()
    n = len(top)
    n_cols = 2
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows), squeeze=False)
    axes_flat = axes.reshape(-1)

    for i, ct in enumerate(top):
        ax = axes_flat[i]
        g = donor_fraction[donor_fraction["cell_type"] == ct]
        x = g["age"].to_numpy(dtype=float)
        y = (g["fraction"].to_numpy(dtype=float)) * 100.0
        ax.scatter(x, y, s=20, alpha=0.75, color=PALETTE["primary"])
        if len(g) >= 2:
            lr = stats.linregress(x, y)
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            ys = lr.intercept + lr.slope * xs
            line_color = PALETTE["secondary"] if lr.slope >= 0 else PALETTE["danger"]
            ax.plot(xs, ys, linewidth=2, color=line_color)
        stat = trends.loc[trends["cell_type"] == ct].iloc[0]
        ax.set_title(
            f"{ct}\nrho={stat['spearman_rho']:.2f}, FDR={stat['spearman_fdr']:.2e}, "
            f"n={int(stat['n_donors'])}"
        )
        ax.set_xlabel("Age")
        ax.set_ylabel("Donor-level fraction (%)")
        style_axis(ax, grid="y")

    for ax in axes_flat[n:]:
        ax.axis("off")

    finalize_and_save(fig, path, dpi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--fig-fractions", required=True)
    ap.add_argument("--fig-trends", required=True)
    ap.add_argument("--table-donor-fractions", required=True)
    ap.add_argument("--table-trends", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ccfg = cfg.get("composition_age", {})
    analysis_cfg = cfg.get("analysis", {})

    dpi = int(ccfg.get("dpi", cfg.get("report", {}).get("dpi", 150)))
    top_n_celltypes = int(ccfg.get("top_n_celltypes", 12))
    top_n_trends = int(ccfg.get("top_n_trends", 8))
    min_cells_per_donor = int(ccfg.get("min_cells_per_donor", 200))
    min_donors_per_celltype = int(ccfg.get("min_donors_per_celltype", 8))
    age_bins = ccfg.get("age_bins", analysis_cfg.get("age_bins", [25, 35, 45, 55, 65, 75, 85]))

    ensure_dir(Path(args.fig_fractions).parent)
    ensure_dir(Path(args.fig_trends).parent)
    ensure_dir(Path(args.table_donor_fractions).parent)
    ensure_dir(Path(args.table_trends).parent)
    apply_publication_style(dpi=dpi)

    obs, cols = _prepare_obs(args.inp, cfg)
    if obs.empty:
        pd.DataFrame(columns=["donor_id", "age", "cell_type", "n_cells", "total_cells", "fraction"]).to_csv(
            args.table_donor_fractions,
            index=False,
        )
        pd.DataFrame(
            columns=[
                "cell_type",
                "n_donors",
                "mean_fraction",
                "spearman_rho",
                "spearman_pvalue",
                "slope_per_year",
                "slope_pvalue",
                "r_squared",
            ]
        ).to_csv(args.table_trends, index=False)
        missing = [k for k, v in cols.items() if v is None]
        msg = f"Missing required obs columns: {', '.join(missing)}"
        _save_placeholder(Path(args.fig_fractions), "Cell-type Composition by Age Bin", msg, dpi)
        _save_placeholder(Path(args.fig_trends), "Top Cell-type Age Trends", msg, dpi)
        return

    donor_fraction = _build_donor_fraction_table(obs, min_cells_per_donor=min_cells_per_donor)
    trends = _compute_trends(donor_fraction, min_donors_per_celltype=min_donors_per_celltype)

    donor_fraction.to_csv(args.table_donor_fractions, index=False)
    trends.to_csv(args.table_trends, index=False)

    _plot_fraction_by_age_bin(
        donor_fraction,
        path=Path(args.fig_fractions),
        age_bins=age_bins,
        top_n_celltypes=top_n_celltypes,
        dpi=dpi,
    )
    _plot_top_trends(
        donor_fraction,
        trends,
        path=Path(args.fig_trends),
        top_n_trends=top_n_trends,
        dpi=dpi,
    )


if __name__ == "__main__":
    main()
