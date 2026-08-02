from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .biological_replicates import attach_biological_replicates
from .plot_style import (
    PALETTE,
    apply_publication_style,
    categorical_palette,
    finalize_and_save,
    save_placeholder,
    style_axis,
)
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


def _bootstrap_rho_slope_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int,
    ci: float,
    seed: int,
) -> dict[str, float]:
    out = {
        "spearman_rho_ci_low": np.nan,
        "spearman_rho_ci_high": np.nan,
        "slope_per_year_ci_low": np.nan,
        "slope_per_year_ci_high": np.nan,
    }
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    n = x.size
    if n < 4 or int(n_boot) <= 1:
        return out

    rng = np.random.default_rng(int(seed))
    rho_vals: list[float] = []
    slope_vals: list[float] = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        xb = x[idx]
        yb = y[idx]
        if np.unique(xb).size < 3 or np.unique(yb).size < 3:
            continue
        try:
            rho, _ = stats.spearmanr(xb, yb)
            slope = stats.linregress(xb, yb).slope
        except Exception:
            continue
        if np.isfinite(rho):
            rho_vals.append(float(rho))
        if np.isfinite(slope):
            slope_vals.append(float(slope))

    alpha = (1.0 - float(ci)) / 2.0
    if len(rho_vals) > 2:
        out["spearman_rho_ci_low"] = float(np.quantile(rho_vals, alpha))
        out["spearman_rho_ci_high"] = float(np.quantile(rho_vals, 1.0 - alpha))
    if len(slope_vals) > 2:
        out["slope_per_year_ci_low"] = float(np.quantile(slope_vals, alpha))
        out["slope_per_year_ci_high"] = float(np.quantile(slope_vals, 1.0 - alpha))
    return out


def _mode_or_na(values: pd.Series) -> str:
    mode = values.mode(dropna=True)
    if mode.empty:
        return "NA"
    return str(mode.iloc[0])


def _resolve_covariate_cols(columns: list[str], cfg: dict) -> list[str]:
    ccfg = cfg.get("composition_age", {})
    requested = ccfg.get("covariate_cols")
    if isinstance(requested, list) and requested:
        candidates = [str(c) for c in requested]
    else:
        candidates = ["sex", "batch"]
    cols = set(columns)
    return [c for c in candidates if c in cols]


def _build_covariate_design(df: pd.DataFrame) -> np.ndarray | None:
    if df.empty:
        return None

    mats: list[np.ndarray] = []
    for col in df.columns:
        s = df[col]
        num = pd.to_numeric(s, errors="coerce")
        numeric_ratio = float(num.notna().mean())
        if numeric_ratio > 0.95:
            if num.notna().any():
                fill = float(num.median(skipna=True))
            else:
                fill = 0.0
            arr = num.fillna(fill).to_numpy(dtype=float).reshape(-1, 1)
            mats.append(arr)
            continue

        cat = s.astype(str).fillna("NA")
        dummies = pd.get_dummies(cat, prefix=str(col), drop_first=True)
        if dummies.shape[1] > 0:
            mats.append(dummies.to_numpy(dtype=float))

    if not mats:
        return None
    x = np.concatenate(mats, axis=1)
    intercept = np.ones((x.shape[0], 1), dtype=float)
    return np.concatenate([intercept, x], axis=1)


def _residualize_vector(y: np.ndarray, x: np.ndarray | None) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if x is None or x.shape[0] != y.shape[0] or x.shape[1] == 0:
        return y - float(np.mean(y))
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    res = y - x @ beta
    if not np.isfinite(res).all() or float(np.std(res)) == 0.0:
        return y - float(np.mean(y))
    return res


def _prepare_xy_for_stats(
    g: pd.DataFrame,
    covariate_cols: list[str],
    adjust_covariates: bool,
) -> tuple[np.ndarray, np.ndarray, bool]:
    x = g["age"].to_numpy(dtype=float)
    y = g["fraction"].to_numpy(dtype=float)
    if not adjust_covariates or not covariate_cols:
        return x, y, False

    design = _build_covariate_design(g[covariate_cols].copy())
    if design is None:
        return x, y, False
    x_res = _residualize_vector(x, design)
    y_res = _residualize_vector(y, design)
    return x_res, y_res, True


def _prepare_obs(inp_h5ad: str, cfg: dict) -> tuple[pd.DataFrame, dict]:
    ccfg = cfg.get("composition_age", {})
    adata = ad.read_h5ad(inp_h5ad, backed="r")
    attach_biological_replicates(adata, cfg)
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
        else _first_existing(
            all_cols, ["donor_id", "donor_id_y", "donor_id_x", "donor", "sample_id"]
        )
    )
    celltype_col = (
        cfg_celltype_col
        if cfg_celltype_col in obs.columns
        else _first_existing(
            all_cols, ["cell_type", "majority_voting", "predicted_labels", "leiden"]
        )
    )

    required = {
        "age_col": age_col,
        "donor_col": donor_col,
        "celltype_col": celltype_col,
    }
    covariate_cols = _resolve_covariate_cols(all_cols, cfg)
    required["covariate_cols"] = covariate_cols
    if age_col is None or donor_col is None or celltype_col is None:
        return pd.DataFrame(), required

    keep_cols = [age_col, donor_col, celltype_col] + [
        c for c in covariate_cols if c not in {age_col, donor_col, celltype_col}
    ]
    out = obs[keep_cols].copy()
    out = out.rename(columns={age_col: "age", donor_col: "donor_id", celltype_col: "cell_type"})
    out["age"] = pd.to_numeric(out["age"], errors="coerce")
    out["donor_id"] = out["donor_id"].astype(str)
    out["cell_type"] = out["cell_type"].astype(str)
    out = out.dropna(subset=["age", "donor_id", "cell_type"])
    return out, required


def _build_donor_fraction_table(
    obs: pd.DataFrame, min_cells_per_donor: int, covariate_cols: list[str]
) -> pd.DataFrame:
    agg_map: dict[str, tuple[str, str | callable]] = {
        "age": ("age", "median"),
        "total_cells": ("cell_type", "size"),
    }
    for c in covariate_cols:
        if c in obs.columns:
            agg_map[c] = (c, _mode_or_na)

    totals = obs.groupby("donor_id", observed=False).agg(**agg_map)
    totals = totals[totals["total_cells"] >= min_cells_per_donor].reset_index()
    if totals.empty:
        return pd.DataFrame(
            columns=["donor_id", "age", "cell_type", "n_cells", "total_cells", "fraction"]
        )

    kept = obs[obs["donor_id"].isin(totals["donor_id"])]
    counts = (
        kept.groupby(["donor_id", "cell_type"], observed=False).size().reset_index(name="n_cells")
    )
    out = counts.merge(totals, on="donor_id", how="left")
    out["fraction"] = out["n_cells"] / out["total_cells"]
    keep = (
        ["donor_id", "age"]
        + [c for c in covariate_cols if c in out.columns]
        + ["cell_type", "n_cells", "total_cells", "fraction"]
    )
    return out[keep]


def _compute_trends(
    donor_fraction: pd.DataFrame,
    min_donors_per_celltype: int,
    covariate_cols: list[str],
    adjust_covariates: bool,
    bootstrap_iterations: int,
    bootstrap_ci: float,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict] = []
    adjusted_for = ",".join(covariate_cols) if (adjust_covariates and covariate_cols) else ""
    for i, (ct, g) in enumerate(donor_fraction.groupby("cell_type", observed=False)):
        if len(g) < min_donors_per_celltype:
            continue
        x, y, used_adjustment = _prepare_xy_for_stats(
            g, covariate_cols=covariate_cols, adjust_covariates=adjust_covariates
        )
        rho, p = stats.spearmanr(x, y)
        lr = stats.linregress(x, y)
        if not np.isfinite(rho) or not np.isfinite(p):
            continue
        ci_stats = _bootstrap_rho_slope_ci(
            x=x,
            y=y,
            n_boot=bootstrap_iterations,
            ci=bootstrap_ci,
            seed=seed + (i * 97),
        )
        rows.append(
            {
                "cell_type": ct,
                "n_donors": int(len(g)),
                "mean_fraction": float(g["fraction"].mean()),
                "spearman_rho": float(rho),
                "spearman_pvalue": float(p),
                "spearman_rho_ci_low": ci_stats["spearman_rho_ci_low"],
                "spearman_rho_ci_high": ci_stats["spearman_rho_ci_high"],
                "slope_per_year": float(lr.slope),
                "slope_per_year_ci_low": ci_stats["slope_per_year_ci_low"],
                "slope_per_year_ci_high": ci_stats["slope_per_year_ci_high"],
                "slope_per_10y": float(lr.slope * 10.0),
                "slope_per_10y_ci_low": float(ci_stats["slope_per_year_ci_low"] * 10.0)
                if np.isfinite(ci_stats["slope_per_year_ci_low"])
                else np.nan,
                "slope_per_10y_ci_high": float(ci_stats["slope_per_year_ci_high"] * 10.0)
                if np.isfinite(ci_stats["slope_per_year_ci_high"])
                else np.nan,
                "slope_pvalue": float(lr.pvalue),
                "r_squared": float(lr.rvalue**2),
                "adjusted": bool(used_adjustment),
                "adjusted_for": adjusted_for if used_adjustment else "",
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
                "spearman_rho_ci_low",
                "spearman_rho_ci_high",
                "spearman_fdr",
                "slope_per_year",
                "slope_per_year_ci_low",
                "slope_per_year_ci_high",
                "slope_per_10y",
                "slope_per_10y_ci_low",
                "slope_per_10y_ci_high",
                "slope_pvalue",
                "slope_fdr",
                "r_squared",
                "adjusted",
                "adjusted_for",
                "direction",
                "fdr_significant",
            ]
        )
    out = pd.DataFrame(rows)
    out["spearman_fdr"] = _bh_fdr(out["spearman_pvalue"].to_numpy())
    out["slope_fdr"] = _bh_fdr(out["slope_pvalue"].to_numpy())
    out["direction"] = np.where(out["slope_per_year"] > 0, "increase_with_age", "decrease_with_age")
    out["fdr_significant"] = out["spearman_fdr"] < 0.05
    return out.sort_values(
        ["spearman_fdr", "spearman_pvalue", "spearman_rho"], ascending=[True, True, False]
    )


def _plot_fraction_by_age_bin(
    donor_fraction: pd.DataFrame,
    path: Path,
    age_bins: list[float],
    top_n_celltypes: int,
    min_donors_per_bin: int,
    drop_sparse_bins: bool,
    sparse_bin_warn_threshold: int,
    sparse_bin_marker: str,
    dpi: int,
) -> None:
    if donor_fraction.empty:
        _save_placeholder(
            path, "Cell-type Composition by Age Bin", "No donor-level fractions available.", dpi
        )
        return

    df = donor_fraction.copy()
    bins = sorted(set(float(v) for v in age_bins))
    if len(bins) < 2:
        _save_placeholder(
            path, "Cell-type Composition by Age Bin", "Need at least two age bin edges.", dpi
        )
        return

    labels = [f"{int(bins[i])}-{int(bins[i + 1])}" for i in range(len(bins) - 1)]
    df["age_bin"] = pd.cut(df["age"], bins=bins, labels=labels, include_lowest=True, right=False)
    df = df.dropna(subset=["age_bin"])
    if df.empty:
        _save_placeholder(
            path,
            "Cell-type Composition by Age Bin",
            "No donors fall inside configured age bins.",
            dpi,
        )
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
    donors_per_bin = (
        df.groupby("age_bin", observed=False)["donor_id"]
        .nunique()
        .reindex(pivot.index)
        .fillna(0)
        .astype(int)
    )
    if min_donors_per_bin > 0 and drop_sparse_bins:
        keep_bins = donors_per_bin >= int(min_donors_per_bin)
        pivot = pivot.loc[keep_bins]
        donors_per_bin = donors_per_bin.loc[keep_bins]
        if pivot.empty:
            _save_placeholder(
                path,
                "Cell-type Composition by Age Bin",
                f"No age bins with donors >= {int(min_donors_per_bin)}.",
                dpi,
            )
            return
    sparse_warn_threshold = int(max(sparse_bin_warn_threshold, 0))
    marker = str(sparse_bin_marker).strip() or "*"
    sparse_mask = (
        donors_per_bin < sparse_warn_threshold
        if sparse_warn_threshold > 0
        else pd.Series(False, index=donors_per_bin.index)
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(pivot.index))
    bottom = np.zeros(len(pivot.index), dtype=float)
    cols = list(pivot.columns)
    bar_colors = categorical_palette(cols, other_label="Other", other_color="#B0B0B0")
    for col, col_color in zip(cols, bar_colors):
        y = (pivot[col] * 100.0).to_numpy(dtype=float)
        ax.bar(x, y, bottom=bottom, label=col, width=0.85, color=col_color)
        bottom += y

    ax.set_xticks(x)
    xt = []
    donors_counts = donors_per_bin.to_numpy()
    for b, n in zip(pivot.index, donors_counts):
        warn = marker if (sparse_warn_threshold > 0 and int(n) < sparse_warn_threshold) else ""
        xt.append(f"{str(b)}{warn}\n(n={int(n)})")
    ax.set_xticklabels(xt, rotation=30, ha="right")
    for tick, n in zip(ax.get_xticklabels(), donors_counts):
        if sparse_warn_threshold > 0 and int(n) < sparse_warn_threshold:
            tick.set_color(PALETTE["muted"])
    ax.set_ylabel("Mean fraction across donors (%)")
    ax.set_xlabel("Age bin")
    ax.set_ylim(0, 100)
    title = "Cell-type composition across age bins\nDonor-level mean fractions"
    if min_donors_per_bin > 0 and drop_sparse_bins:
        title += f" (bins with n>={int(min_donors_per_bin)})"
    elif sparse_warn_threshold > 0 and bool(np.any(sparse_mask.to_numpy(dtype=bool))):
        title += f" ({marker} low support: n<{sparse_warn_threshold})"
    ax.set_title(title)
    style_axis(ax, grid="y")
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    finalize_and_save(fig, path, dpi)


def _plot_top_trends(
    donor_fraction: pd.DataFrame,
    trends: pd.DataFrame,
    path: Path,
    top_n_trends: int,
    min_donors_for_panel: int,
    significant_only: bool,
    dpi: int,
) -> None:
    if donor_fraction.empty or trends.empty:
        _save_placeholder(path, "Top Cell-type Age Trends", "No trend statistics available.", dpi)
        return

    candidates = trends.copy()
    if significant_only:
        candidates = candidates[candidates["fdr_significant"]].copy()
    if min_donors_for_panel > 0:
        candidates = candidates[candidates["n_donors"] >= int(min_donors_for_panel)].copy()

    if candidates.empty:
        filters = []
        if significant_only:
            filters.append("FDR<0.05")
        if min_donors_for_panel > 0:
            filters.append(f"n_donors>={min_donors_for_panel}")
        suffix = f" ({', '.join(filters)})" if filters else ""
        _save_placeholder(
            path, "Top Cell-type Age Trends", f"No trends satisfy panel filters{suffix}.", dpi
        )
        return

    top = candidates.head(top_n_trends)["cell_type"].tolist()
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
        adjusted_note = ""
        if bool(stat.get("adjusted", False)):
            adjusted_for = str(stat.get("adjusted_for", "")).strip()
            adjusted_note = f", adj={adjusted_for}" if adjusted_for else ", adj"
        ax.set_title(
            f"{ct}\nrho={stat['spearman_rho']:.2f}, FDR={stat['spearman_fdr']:.2e}, "
            f"n={int(stat['n_donors'])}{adjusted_note}"
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
    min_donors_for_panel = int(ccfg.get("min_donors_for_panel", max(20, min_donors_per_celltype)))
    plot_only_fdr_significant = bool(ccfg.get("plot_only_fdr_significant", True))
    adjust_covariates = bool(ccfg.get("adjust_covariates", True))
    min_donors_per_age_bin = int(ccfg.get("min_donors_per_age_bin", 12))
    drop_sparse_age_bins = bool(ccfg.get("drop_sparse_age_bins", True))
    sparse_age_bin_warn_threshold = int(
        ccfg.get("sparse_age_bin_warn_threshold", min_donors_per_age_bin)
    )
    sparse_age_bin_marker = str(ccfg.get("sparse_age_bin_marker", "*"))
    bootstrap_iterations = int(
        ccfg.get(
            "bootstrap_iterations", cfg.get("age_prediction", {}).get("bootstrap_iterations", 2000)
        )
    )
    bootstrap_ci = float(
        ccfg.get("bootstrap_ci", cfg.get("age_prediction", {}).get("bootstrap_ci", 0.95))
    )
    seed = int(cfg.get("run", {}).get("seed", 42))
    age_bins = ccfg.get("age_bins", analysis_cfg.get("age_bins", [25, 35, 45, 55, 65, 75, 85]))
    bootstrap_iterations = max(0, bootstrap_iterations)
    bootstrap_ci = float(np.clip(bootstrap_ci, 0.5, 0.999))

    ensure_dir(Path(args.fig_fractions).parent)
    ensure_dir(Path(args.fig_trends).parent)
    ensure_dir(Path(args.table_donor_fractions).parent)
    ensure_dir(Path(args.table_trends).parent)
    apply_publication_style(dpi=dpi)

    obs, cols = _prepare_obs(args.inp, cfg)
    print(
        f"[composition_age] prepared {len(obs):,} cells for aggregation",
        flush=True,
    )
    if obs.empty:
        pd.DataFrame(
            columns=["donor_id", "age", "cell_type", "n_cells", "total_cells", "fraction"]
        ).to_csv(
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
                "spearman_rho_ci_low",
                "spearman_rho_ci_high",
                "slope_per_year",
                "slope_per_year_ci_low",
                "slope_per_year_ci_high",
                "slope_per_10y",
                "slope_per_10y_ci_low",
                "slope_per_10y_ci_high",
                "slope_pvalue",
                "r_squared",
            ]
        ).to_csv(args.table_trends, index=False)
        missing = [k for k, v in cols.items() if v is None]
        msg = f"Missing required obs columns: {', '.join(missing)}"
        _save_placeholder(Path(args.fig_fractions), "Cell-type Composition by Age Bin", msg, dpi)
        _save_placeholder(Path(args.fig_trends), "Top Cell-type Age Trends", msg, dpi)
        return

    covariate_cols = list(cols.get("covariate_cols", []))
    donor_fraction = _build_donor_fraction_table(
        obs,
        min_cells_per_donor=min_cells_per_donor,
        covariate_cols=covariate_cols,
    )
    trends = _compute_trends(
        donor_fraction,
        min_donors_per_celltype=min_donors_per_celltype,
        covariate_cols=covariate_cols,
        adjust_covariates=adjust_covariates,
        bootstrap_iterations=bootstrap_iterations,
        bootstrap_ci=bootstrap_ci,
        seed=seed,
    )
    print(
        "[composition_age] "
        f"computed {len(trends):,} cell-type trend tests from "
        f"{donor_fraction['donor_id'].nunique():,} replicates",
        flush=True,
    )

    grouping_col = str(ccfg.get("donor_col", "donor_id"))
    donor_fraction_out = donor_fraction.rename(columns={"donor_id": grouping_col})
    trends["grouping_id_column"] = grouping_col
    donor_fraction_out.to_csv(args.table_donor_fractions, index=False)
    trends.to_csv(args.table_trends, index=False)

    _plot_fraction_by_age_bin(
        donor_fraction,
        path=Path(args.fig_fractions),
        age_bins=age_bins,
        top_n_celltypes=top_n_celltypes,
        min_donors_per_bin=min_donors_per_age_bin,
        drop_sparse_bins=drop_sparse_age_bins,
        sparse_bin_warn_threshold=sparse_age_bin_warn_threshold,
        sparse_bin_marker=sparse_age_bin_marker,
        dpi=dpi,
    )
    _plot_top_trends(
        donor_fraction,
        trends,
        path=Path(args.fig_trends),
        top_n_trends=top_n_trends,
        min_donors_for_panel=min_donors_for_panel,
        significant_only=plot_only_fdr_significant,
        dpi=dpi,
    )


if __name__ == "__main__":
    main()
