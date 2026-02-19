from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .plot_style import PALETTE, apply_publication_style, finalize_and_save, save_placeholder, style_axis
from .utils import ensure_dir, load_config

SUMMARY_COLUMNS = [
    "reference_model",
    "reference_params",
    "candidate_model",
    "candidate_params",
    "candidate_is_best",
    "candidate_is_best_ridge",
    "n_folds_compared",
    "winner_count_candidate",
    "winner_count_reference",
    "tie_count",
    "win_fraction_candidate",
    "reference_mae_mean",
    "candidate_mae_mean",
    "delta_mae_mean",
    "delta_mae_median",
    "delta_mae_boot_ci_low",
    "delta_mae_boot_ci_high",
    "delta_mae_boot_prob_gt0",
]


def _save_placeholder(path: Path, title: str, message: str, dpi: int) -> None:
    save_placeholder(path=path, title=title, message=message, dpi=dpi)


def _first_existing(columns: list[str], candidates: list[str]) -> str | None:
    cols = set(columns)
    for key in candidates:
        if key in cols:
            return key
    return None


def _write_empty_outputs(
    fig_pred: Path,
    fig_mae: Path,
    table_pred: Path,
    table_metrics: Path,
    table_summary: Path,
    dpi: int,
    msg: str,
) -> None:
    _save_placeholder(fig_pred, "Age Prediction: Observed vs Predicted", msg, dpi)
    _save_placeholder(fig_mae, "Age Prediction: MAE by Cell Type", msg, dpi)
    pd.DataFrame(
        columns=[
            "donor_id",
            "cell_type",
            "age_true",
            "age_pred",
            "residual",
            "abs_error",
            "within_5y",
            "n_cells",
            "fold",
        ]
    ).to_csv(table_pred, index=False)
    pd.DataFrame(
        columns=[
            "model",
            "family",
            "params",
            "alpha",
            "mae",
            "median_ae",
            "rmse",
            "r2",
            "pearson_r",
            "spearman_rho",
            "within_5y_fraction",
            "baseline_mae",
            "delta_mae_vs_baseline",
            "delta_mae_vs_best_ridge",
            "n_samples",
            "n_donors",
            "n_celltypes",
            "n_splits",
            "is_best_ridge",
            "is_best",
        ]
    ).to_csv(table_metrics, index=False)
    pd.DataFrame(columns=SUMMARY_COLUMNS).to_csv(table_summary, index=False)


def _prepare_aggregated_table(inp_h5ad: str, cfg: dict) -> tuple[pd.DataFrame, list[str], str]:
    pcfg = cfg.get("age_prediction", {})
    latent_key = str(pcfg.get("latent_key", "X_scVI"))
    max_cells = int(pcfg.get("max_cells_for_aggregation", 300000))
    min_cells_per_group = int(pcfg.get("min_cells_per_group", 50))
    seed = int(cfg.get("run", {}).get("seed", 42))

    adata = ad.read_h5ad(inp_h5ad, backed="r")
    try:
        obs_cols = list(adata.obs.columns)

        age_col = pcfg.get("age_col")
        donor_col = pcfg.get("donor_col")
        celltype_col = pcfg.get("celltype_col")

        if age_col not in obs_cols:
            age_col = _first_existing(obs_cols, ["age"])
        if donor_col not in obs_cols:
            donor_col = _first_existing(obs_cols, ["donor_id", "donor_id_y", "donor_id_x", "donor", "sample_id"])
        if celltype_col not in obs_cols:
            celltype_col = _first_existing(obs_cols, ["cell_type", "majority_voting", "predicted_labels", "leiden"])

        if age_col is None or donor_col is None or celltype_col is None:
            return pd.DataFrame(), [], "Missing required obs columns for age prediction."

        if latent_key not in adata.obsm:
            return pd.DataFrame(), [], f"Missing latent representation in obsm: {latent_key}"

        obs = adata.obs[[age_col, donor_col, celltype_col]].copy()
        obs.columns = ["age", "donor_id", "cell_type"]
        obs["age"] = pd.to_numeric(obs["age"], errors="coerce")

        valid_mask = obs["age"].notna() & obs["donor_id"].notna() & obs["cell_type"].notna()
        valid_idx = np.flatnonzero(valid_mask.to_numpy())
        if len(valid_idx) == 0:
            return pd.DataFrame(), [], "No valid cells after age/donor/cell_type filtering."

        if len(valid_idx) > max_cells:
            rng = np.random.default_rng(seed)
            valid_idx = np.sort(rng.choice(valid_idx, size=max_cells, replace=False))

        # Prefer indexed extraction to avoid materializing all cells in memory.
        latent_obj = adata.obsm[latent_key]
        try:
            latent = np.asarray(latent_obj[valid_idx], dtype=np.float32)
        except Exception:
            latent = np.asarray(latent_obj, dtype=np.float32)[valid_idx]

        latent_cols = [f"latent_{i}" for i in range(latent.shape[1])]
        latent_df = pd.DataFrame(latent, columns=latent_cols)
        obs_sub = obs.iloc[valid_idx].reset_index(drop=True)
        cell_df = pd.concat([obs_sub, latent_df], axis=1)

        agg_map = {"age": ("age", "median"), "n_cells": ("age", "size")}
        for c in latent_cols:
            agg_map[c] = (c, "mean")

        agg = (
            cell_df.groupby(["donor_id", "cell_type"], observed=False)
            .agg(**agg_map)
            .reset_index()
        )
        agg = agg[agg["n_cells"] >= min_cells_per_group].copy()

        if agg.empty:
            return pd.DataFrame(), [], "No donor-celltype groups pass min_cells_per_group."

        return agg, latent_cols, ""
    finally:
        if hasattr(adata, "file") and getattr(adata, "file", None) is not None:
            adata.file.close()


def _load_model_builders(cfg: dict, latent_cols: list[str], seed: int):
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler

    pcfg = cfg.get("age_prediction", {})
    model_order = [str(m).lower() for m in pcfg.get("model_order", ["ridge", "xgboost"])]
    ridge_alphas = pcfg.get("ridge_alphas", pcfg.get("alphas", [0.1, 1.0, 10.0, 100.0]))
    ridge_alphas = [float(a) for a in ridge_alphas]
    n_jobs = int(pcfg.get("n_jobs", -1))

    try:
        ohe = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        ohe = OneHotEncoder(handle_unknown="ignore", sparse=False)

    preproc_linear = ColumnTransformer(
        transformers=[
            ("num", StandardScaler(), latent_cols),
            ("cat", ohe, ["cell_type"]),
        ],
        remainder="drop",
    )
    preproc_tree = ColumnTransformer(
        transformers=[
            ("num", "passthrough", latent_cols),
            ("cat", ohe, ["cell_type"]),
        ],
        remainder="drop",
    )

    candidates: list[dict] = []
    if "ridge" in model_order:
        for a in ridge_alphas:
            label = f"ridge_alpha_{a:g}"
            candidates.append(
                {
                    "model": "ridge",
                    "family": "linear",
                    "label": label,
                    "params": {"alpha": float(a)},
                    "build": lambda alpha=float(a): Pipeline(
                        steps=[("preproc", preproc_linear), ("reg", Ridge(alpha=alpha))]
                    ),
                }
            )

    if "hist_gbr" in model_order:
        hcfg = pcfg.get("hist_gbr_params", [{"max_depth": 4, "learning_rate": 0.05, "max_iter": 300}])
        for i, p in enumerate(hcfg):
            params = dict(p)
            params.setdefault("random_state", seed)
            label = f"hist_gbr_{i+1}"
            candidates.append(
                {
                    "model": "hist_gbr",
                    "family": "tree",
                    "label": label,
                    "params": params,
                    "build": lambda params=params: Pipeline(
                        steps=[("preproc", preproc_tree), ("reg", HistGradientBoostingRegressor(**params))]
                    ),
                }
            )

    if "xgboost" in model_order:
        try:
            from xgboost import XGBRegressor

            xcfg = pcfg.get(
                "xgboost_params",
                [
                    {
                        "n_estimators": 300,
                        "max_depth": 3,
                        "learning_rate": 0.05,
                        "subsample": 0.8,
                        "colsample_bytree": 0.8,
                        "reg_lambda": 1.0,
                    },
                    {
                        "n_estimators": 500,
                        "max_depth": 4,
                        "learning_rate": 0.05,
                        "subsample": 0.8,
                        "colsample_bytree": 0.8,
                        "reg_lambda": 1.0,
                    },
                ],
            )
            for i, p in enumerate(xcfg):
                params = dict(p)
                params.setdefault("objective", "reg:squarederror")
                params.setdefault("random_state", seed)
                params.setdefault("n_jobs", n_jobs)
                params.setdefault("verbosity", 0)
                label = f"xgboost_{i+1}"
                candidates.append(
                    {
                        "model": "xgboost",
                        "family": "tree",
                        "label": label,
                        "params": params,
                        "build": lambda params=params: Pipeline(
                            steps=[("preproc", preproc_tree), ("reg", XGBRegressor(**params))]
                        ),
                    }
                )
        except ImportError:
            pass

    if "lightgbm" in model_order:
        try:
            from lightgbm import LGBMRegressor

            lcfg = pcfg.get(
                "lightgbm_params",
                [
                    {
                        "n_estimators": 300,
                        "max_depth": -1,
                        "learning_rate": 0.05,
                        "subsample": 0.8,
                        "colsample_bytree": 0.8,
                        "reg_lambda": 1.0,
                    }
                ],
            )
            for i, p in enumerate(lcfg):
                params = dict(p)
                params.setdefault("random_state", seed)
                params.setdefault("n_jobs", n_jobs)
                label = f"lightgbm_{i+1}"
                candidates.append(
                    {
                        "model": "lightgbm",
                        "family": "tree",
                        "label": label,
                        "params": params,
                        "build": lambda params=params: Pipeline(
                            steps=[("preproc", preproc_tree), ("reg", LGBMRegressor(**params))]
                        ),
                    }
                )
        except ImportError:
            pass

    return candidates


def _evaluate_candidates(agg: pd.DataFrame, latent_cols: list[str], cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.model_selection import GroupKFold

    pcfg = cfg.get("age_prediction", {})
    cv_folds = int(pcfg.get("cv_folds", 5))
    seed = int(cfg.get("run", {}).get("seed", 42))

    y = agg["age"].to_numpy(dtype=float)
    groups = agg["donor_id"].astype(str).to_numpy()
    X = agg[latent_cols + ["cell_type"]].copy()

    unique_groups = np.unique(groups)
    if unique_groups.size < 3:
        raise ValueError("Need at least 3 donors for cross-donor CV.")
    n_splits = min(cv_folds, int(unique_groups.size))
    if n_splits < 2:
        raise ValueError("Not enough donor groups for CV.")

    splitter = GroupKFold(n_splits=n_splits)
    splits = list(splitter.split(X, y, groups=groups))

    candidates = _load_model_builders(cfg, latent_cols=latent_cols, seed=seed)
    if not candidates:
        raise ValueError("No valid models available (check model_order and installed libraries).")

    metrics_rows = []
    pred_rows = []
    for cand in candidates:
        model = cand["build"]()
        y_pred = np.full_like(y, np.nan, dtype=float)
        y_base = np.full_like(y, np.nan, dtype=float)
        fold_id = np.full(y.shape[0], -1, dtype=int)

        for fold, (tr, te) in enumerate(splits):
            model.fit(X.iloc[tr], y[tr])
            y_pred[te] = model.predict(X.iloc[te])
            y_base[te] = np.mean(y[tr])
            fold_id[te] = fold

        abs_err = np.abs(y_pred - y)
        mae = float(mean_absolute_error(y, y_pred))
        rmse = float(np.sqrt(mean_squared_error(y, y_pred)))
        r2 = float(r2_score(y, y_pred))
        pearson_r = float(stats.pearsonr(y, y_pred)[0]) if len(y) > 1 else np.nan
        spearman_rho = float(stats.spearmanr(y, y_pred)[0]) if len(y) > 1 else np.nan
        baseline_mae = float(mean_absolute_error(y, y_base))

        metrics_rows.append(
            {
                "model": cand["model"],
                "family": cand["family"],
                "params": json.dumps(cand["params"], sort_keys=True),
                "alpha": float(cand["params"].get("alpha", np.nan)),
                "mae": mae,
                "median_ae": float(np.median(abs_err)),
                "rmse": rmse,
                "r2": r2,
                "pearson_r": pearson_r,
                "spearman_rho": spearman_rho,
                "within_5y_fraction": float(np.mean(abs_err <= 5.0)),
                "baseline_mae": baseline_mae,
                "delta_mae_vs_baseline": baseline_mae - mae,
                "n_samples": int(len(y)),
                "n_donors": int(unique_groups.size),
                "n_celltypes": int(agg["cell_type"].nunique()),
                "n_splits": int(n_splits),
            }
        )

        pred = agg[["donor_id", "cell_type", "age", "n_cells"]].copy()
        pred = pred.rename(columns={"age": "age_true"})
        pred["model"] = cand["model"]
        pred["params"] = json.dumps(cand["params"], sort_keys=True)
        pred["age_pred"] = y_pred
        pred["residual"] = pred["age_pred"] - pred["age_true"]
        pred["abs_error"] = np.abs(pred["residual"])
        pred["within_5y"] = pred["abs_error"] <= 5.0
        pred["fold"] = fold_id
        pred_rows.append(pred)

    metrics = pd.DataFrame(metrics_rows).sort_values("mae", ascending=True).reset_index(drop=True)
    metrics["is_best"] = False
    if not metrics.empty:
        metrics.loc[0, "is_best"] = True

    ridge_mask = metrics["model"].eq("ridge")
    metrics["is_best_ridge"] = False
    metrics["delta_mae_vs_best_ridge"] = np.nan
    if ridge_mask.any():
        ridge_idx = metrics.loc[ridge_mask, "mae"].idxmin()
        ridge_best_mae = float(metrics.loc[ridge_idx, "mae"])
        metrics.loc[ridge_idx, "is_best_ridge"] = True
        metrics["delta_mae_vs_best_ridge"] = ridge_best_mae - metrics["mae"]

    pred_all = pd.concat(pred_rows, axis=0, ignore_index=True)
    best_row = metrics.loc[metrics["is_best"]].iloc[0]
    best_mask = (pred_all["model"] == best_row["model"]) & (pred_all["params"] == best_row["params"])
    pred_all["is_best_model"] = best_mask
    pred_all = pred_all.sort_values(["is_best_model", "model", "donor_id", "cell_type"], ascending=[False, True, True, True])

    return pred_all, metrics


def _bootstrap_mean_ci(values: np.ndarray, n_boot: int, ci: float, seed: int) -> tuple[float, float, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan, np.nan, np.nan
    if vals.size == 1:
        v = float(vals[0])
        return v, v, float(v > 0.0)

    rng = np.random.default_rng(seed)
    idx = rng.integers(0, vals.size, size=(n_boot, vals.size))
    boot_means = vals[idx].mean(axis=1)
    alpha = (1.0 - ci) / 2.0
    low = float(np.quantile(boot_means, alpha))
    high = float(np.quantile(boot_means, 1.0 - alpha))
    prob_gt0 = float(np.mean(boot_means > 0.0))
    return low, high, prob_gt0


def _build_model_comparison_summary(pred_all: pd.DataFrame, metrics: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    pcfg = cfg.get("age_prediction", {})
    n_boot = int(pcfg.get("bootstrap_iterations", 2000))
    ci = float(pcfg.get("bootstrap_ci", 0.95))
    seed = int(cfg.get("run", {}).get("seed", 42))

    if pred_all.empty or metrics.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    pred_all = pred_all.copy()
    pred_all["model_key"] = pred_all["model"] + "||" + pred_all["params"].astype(str)

    fold_mae = (
        pred_all.groupby(["model", "params", "model_key", "fold"], observed=False)["abs_error"]
        .mean()
        .reset_index(name="fold_mae")
    )
    pivot = fold_mae.pivot(index="fold", columns="model_key", values="fold_mae")

    ridge_rows = metrics[metrics["is_best_ridge"]]
    if ridge_rows.empty:
        ref_row = metrics.iloc[0]
    else:
        ref_row = ridge_rows.iloc[0]

    ref_key = f"{ref_row['model']}||{ref_row['params']}"
    if ref_key not in pivot.columns:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    summary_rows = []
    for row in metrics.itertuples(index=False):
        cand_key = f"{row.model}||{row.params}"
        if cand_key not in pivot.columns:
            continue

        pair = pivot[[ref_key, cand_key]].dropna()
        if pair.empty:
            continue
        ref = pair[ref_key].to_numpy(dtype=float)
        cand = pair[cand_key].to_numpy(dtype=float)
        delta = ref - cand  # >0 means candidate better (lower MAE)

        winner_candidate = int(np.sum(delta > 0))
        winner_reference = int(np.sum(delta < 0))
        ties = int(np.sum(np.isclose(delta, 0.0)))
        ci_low, ci_high, prob_gt0 = _bootstrap_mean_ci(delta, n_boot=n_boot, ci=ci, seed=seed + len(summary_rows))

        summary_rows.append(
            {
                "reference_model": ref_row["model"],
                "reference_params": ref_row["params"],
                "candidate_model": row.model,
                "candidate_params": row.params,
                "candidate_is_best": bool(row.is_best),
                "candidate_is_best_ridge": bool(row.is_best_ridge),
                "n_folds_compared": int(pair.shape[0]),
                "winner_count_candidate": winner_candidate,
                "winner_count_reference": winner_reference,
                "tie_count": ties,
                "win_fraction_candidate": float(winner_candidate / float(max(pair.shape[0], 1))),
                "reference_mae_mean": float(np.mean(ref)),
                "candidate_mae_mean": float(np.mean(cand)),
                "delta_mae_mean": float(np.mean(delta)),
                "delta_mae_median": float(np.median(delta)),
                "delta_mae_boot_ci_low": ci_low,
                "delta_mae_boot_ci_high": ci_high,
                "delta_mae_boot_prob_gt0": prob_gt0,
            }
        )

    out = pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
    if out.empty:
        return out
    return out.sort_values(["candidate_is_best", "delta_mae_mean"], ascending=[False, False]).reset_index(drop=True)


def _plot_observed_vs_pred(pred_best: pd.DataFrame, metrics: pd.DataFrame, path: Path, dpi: int) -> None:
    best = metrics.loc[metrics["is_best"]].iloc[0]
    x = pred_best["age_true"].to_numpy(dtype=float)
    y = pred_best["age_pred"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    hb = ax.hexbin(x, y, gridsize=32, cmap="Blues", mincnt=1)
    cbar = fig.colorbar(hb, ax=ax, shrink=0.8)
    cbar.set_label("Donor-celltype groups")
    ax.scatter(x, y, s=10, alpha=0.25, color=PALETTE["primary"])
    mn = float(np.nanmin([x.min(), y.min()]))
    mx = float(np.nanmax([x.max(), y.max()]))
    ax.plot([mn, mx], [mn, mx], linestyle="--", linewidth=2, color=PALETTE["muted"], label="Ideal")
    if len(x) >= 2:
        lr = stats.linregress(x, y)
        xs = np.linspace(mn, mx, 100)
        ys = lr.intercept + lr.slope * xs
        ax.plot(xs, ys, linewidth=2, color=PALETTE["danger"], label="Fit")
    ax.set_xlabel("Observed age")
    ax.set_ylabel("Predicted age (cross-validated)")
    ax.set_title(f"Chronological age prediction (best model: {best['model']})")
    top_models = metrics.head(min(3, metrics.shape[0]))[["model", "mae"]]
    top_txt = "; ".join([f"{r.model}:{r.mae:.2f}" for r in top_models.itertuples(index=False)])
    txt = (
        f"Model={best['model']}\n"
        f"MAE={best['mae']:.2f}, median AE={best['median_ae']:.2f}, RMSE={best['rmse']:.2f}\n"
        f"R2={best['r2']:.3f}, Pearson r={best['pearson_r']:.3f}, within 5y={best['within_5y_fraction']:.2%}\n"
        f"Top MAE: {top_txt}"
    )
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=10)
    style_axis(ax, grid="y")
    ax.legend(frameon=False, loc="lower right")
    finalize_and_save(fig, path, dpi)


def _plot_mae_by_celltype(pred_best: pd.DataFrame, best_row: pd.Series, path: Path, top_n: int, dpi: int) -> None:
    ct = (
        pred_best.groupby("cell_type", observed=False)
        .agg(
            mae=("abs_error", "mean"),
            median_ae=("abs_error", "median"),
            n_groups=("abs_error", "size"),
            within_5y_fraction=("within_5y", "mean"),
        )
        .sort_values("mae", ascending=False)
    )
    if ct.empty:
        _save_placeholder(path, "Age Prediction: MAE by Cell Type", "No prediction rows available.", dpi)
        return

    ct = ct.head(top_n).sort_values("mae", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.barh(ct.index, ct["mae"].to_numpy(), color=PALETTE["primary"])
    ax.set_xlabel("Mean absolute error (years)")
    ax.set_ylabel("Cell type")
    ax.set_title(f"Cross-validated age prediction error by cell type\nBest model: {best_row['model']}")
    for i, (_, row) in enumerate(ct.iterrows()):
        ax.text(
            row["mae"],
            i,
            f"  n={int(row['n_groups'])}, medAE={row['median_ae']:.1f}, within5y={row['within_5y_fraction']:.1%}",
            va="center",
            fontsize=8,
        )
    ax.bar_label(bars, labels=[f"{v:.1f}" for v in ct["mae"].to_numpy()], padding=3, fontsize=8)
    style_axis(ax, grid="x")
    finalize_and_save(fig, path, dpi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--fig-pred", required=True)
    ap.add_argument("--fig-mae", required=True)
    ap.add_argument("--table-pred", required=True)
    ap.add_argument("--table-metrics", required=True)
    ap.add_argument("--table-summary", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    pcfg = cfg.get("age_prediction", {})
    dpi = int(pcfg.get("dpi", cfg.get("report", {}).get("dpi", 150)))
    top_n_celltypes = int(pcfg.get("top_n_celltypes_plot", 12))

    fig_pred = Path(args.fig_pred)
    fig_mae = Path(args.fig_mae)
    table_pred = Path(args.table_pred)
    table_metrics = Path(args.table_metrics)
    table_summary = Path(args.table_summary)

    ensure_dir(fig_pred.parent)
    ensure_dir(fig_mae.parent)
    ensure_dir(table_pred.parent)
    ensure_dir(table_metrics.parent)
    ensure_dir(table_summary.parent)
    apply_publication_style(dpi=dpi)

    try:
        import sklearn  # noqa: F401
    except ImportError:
        _write_empty_outputs(
            fig_pred=fig_pred,
            fig_mae=fig_mae,
            table_pred=table_pred,
            table_metrics=table_metrics,
            table_summary=table_summary,
            dpi=dpi,
            msg="scikit-learn is not installed in the current environment.",
        )
        return

    agg, latent_cols, err = _prepare_aggregated_table(args.inp, cfg)
    if agg.empty:
        _write_empty_outputs(
            fig_pred=fig_pred,
            fig_mae=fig_mae,
            table_pred=table_pred,
            table_metrics=table_metrics,
            table_summary=table_summary,
            dpi=dpi,
            msg=err or "No aggregated donor-celltype rows available.",
        )
        return

    try:
        pred_all, metrics = _evaluate_candidates(agg=agg, latent_cols=latent_cols, cfg=cfg)
    except Exception as exc:
        _write_empty_outputs(
            fig_pred=fig_pred,
            fig_mae=fig_mae,
            table_pred=table_pred,
            table_metrics=table_metrics,
            table_summary=table_summary,
            dpi=dpi,
            msg=f"Age prediction failed: {exc}",
        )
        return

    summary = _build_model_comparison_summary(pred_all=pred_all, metrics=metrics, cfg=cfg)
    pred_all.to_csv(table_pred, index=False)
    metrics.to_csv(table_metrics, index=False)
    summary.to_csv(table_summary, index=False)
    best_row = metrics.loc[metrics["is_best"]].iloc[0]
    pred_best = pred_all[pred_all["is_best_model"]].copy()
    _plot_observed_vs_pred(pred_best, metrics, fig_pred, dpi)
    _plot_mae_by_celltype(pred_best, best_row=best_row, path=fig_mae, top_n=top_n_celltypes, dpi=dpi)


if __name__ == "__main__":
    main()
