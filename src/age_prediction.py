from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .annotation_mapping import attach_analysis_cell_types
from .biological_replicates import attach_biological_replicates
from .plot_style import (
    PALETTE,
    apply_publication_style,
    finalize_and_save,
    save_placeholder,
    style_axis,
)
from .scientific_guardrails import require_placeholder_permission
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


def _weighted_capped_allocation(weights: np.ndarray, caps: np.ndarray, total: int) -> np.ndarray:
    caps = np.asarray(caps, dtype=int)
    out = np.zeros_like(caps, dtype=int)
    if caps.size == 0:
        return out

    total = int(max(total, 0))
    if total == 0:
        return out
    budget = min(total, int(caps.sum()))
    if budget <= 0:
        return out

    w = np.asarray(weights, dtype=float)
    w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
    if float(w.sum()) <= 0:
        w = caps.astype(float)
    if float(w.sum()) <= 0:
        w = (caps > 0).astype(float)

    quota = (w / float(w.sum())) * float(budget)
    add = np.floor(quota).astype(int)
    add = np.minimum(add, caps)
    out += add
    remaining = budget - int(out.sum())
    if remaining <= 0:
        return out

    frac = quota - np.floor(quota)
    order = np.argsort(-frac)
    for idx in order:
        if remaining <= 0:
            break
        if out[idx] >= caps[idx]:
            continue
        out[idx] += 1
        remaining -= 1

    if remaining <= 0:
        return out

    cap_left = caps - out
    order = np.argsort(-cap_left)
    for idx in order:
        if remaining <= 0:
            break
        if cap_left[idx] <= 0:
            continue
        take = min(int(cap_left[idx]), remaining)
        out[idx] += take
        remaining -= take
    return out


def _allocate_group_samples(
    group_sizes: pd.Series, max_cells: int, min_cells_per_group: int
) -> pd.Series:
    alloc = pd.Series(0, index=group_sizes.index, dtype=int)
    if group_sizes.empty or max_cells <= 0:
        return alloc

    min_cells_per_group = max(int(min_cells_per_group), 1)
    large_mask = group_sizes >= min_cells_per_group
    mandatory = pd.Series(0, index=group_sizes.index, dtype=int)
    mandatory.loc[large_mask] = min_cells_per_group
    mandatory_total = int(mandatory.sum())

    if mandatory_total >= max_cells:
        weights = group_sizes.where(large_mask, 0).to_numpy(dtype=float)
        caps = group_sizes.where(large_mask, 0).to_numpy(dtype=int)
        alloc_vals = _weighted_capped_allocation(weights=weights, caps=caps, total=max_cells)
        return pd.Series(alloc_vals, index=group_sizes.index, dtype=int)

    alloc += mandatory
    remaining_budget = max_cells - int(alloc.sum())
    if remaining_budget <= 0:
        return alloc

    caps_left = (group_sizes - alloc).clip(lower=0).to_numpy(dtype=int)
    weights = caps_left.astype(float)
    extra = _weighted_capped_allocation(weights=weights, caps=caps_left, total=remaining_budget)
    alloc += pd.Series(extra, index=group_sizes.index, dtype=int)
    return alloc


def _stratified_sample_indices(
    valid_idx: np.ndarray,
    donor_values: np.ndarray,
    celltype_values: np.ndarray,
    max_cells: int,
    min_cells_per_group: int,
    seed: int,
) -> np.ndarray:
    if len(valid_idx) <= max_cells:
        out = np.asarray(valid_idx, dtype=int)
        out.sort()
        return out

    df = pd.DataFrame(
        {
            "obs_idx": np.asarray(valid_idx, dtype=int),
            "donor_id": np.asarray(donor_values, dtype=str),
            "cell_type": np.asarray(celltype_values, dtype=str),
        }
    )
    group_sizes = df.groupby(["donor_id", "cell_type"], observed=False).size()
    alloc = _allocate_group_samples(
        group_sizes=group_sizes,
        max_cells=int(max_cells),
        min_cells_per_group=int(min_cells_per_group),
    )

    rng = np.random.default_rng(seed)
    chosen: list[np.ndarray] = []
    grouped = df.groupby(["donor_id", "cell_type"], observed=False)["obs_idx"]
    for key, obs_idx_series in grouped:
        take = int(alloc.loc[key]) if key in alloc.index else 0
        if take <= 0:
            continue
        arr = obs_idx_series.to_numpy(dtype=int)
        if take >= arr.size:
            chosen.append(arr)
        else:
            chosen.append(rng.choice(arr, size=take, replace=False))

    if not chosen:
        fallback = np.asarray(valid_idx, dtype=int)
        fallback = rng.choice(fallback, size=max_cells, replace=False)
        fallback.sort()
        return fallback

    out = np.concatenate(chosen).astype(int, copy=False)
    if out.size > max_cells:
        out = rng.choice(out, size=max_cells, replace=False)
    out.sort()
    return out


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
            "evaluation_level",
            "donor_id",
            "cell_type",
            "age_true",
            "age_pred",
            "age_base",
            "residual",
            "abs_error",
            "within_5y",
            "n_cells",
            "n_groups",
            "fold",
            "model",
            "params",
            "selected_model",
            "selected_params",
            "is_best_model",
        ]
    ).to_csv(table_pred, index=False)
    pd.DataFrame(
        columns=[
            "model",
            "family",
            "params",
            "alpha",
            "mae",
            "mae_ci_low",
            "mae_ci_high",
            "median_ae",
            "rmse",
            "r2",
            "r2_ci_low",
            "r2_ci_high",
            "pearson_r",
            "spearman_rho",
            "within_5y_fraction",
            "baseline_mae",
            "delta_mae_vs_baseline",
            "delta_mae_vs_best_ridge",
            "mae_donor_level",
            "mae_donor_ci_low",
            "mae_donor_ci_high",
            "median_ae_donor_level",
            "rmse_donor_level",
            "r2_donor_level",
            "r2_donor_ci_low",
            "r2_donor_ci_high",
            "pearson_r_donor_level",
            "spearman_rho_donor_level",
            "within_5y_fraction_donor_level",
            "baseline_mae_donor_level",
            "delta_mae_vs_baseline_donor_level",
            "mae_group_level",
            "median_ae_group_level",
            "rmse_group_level",
            "r2_group_level",
            "pearson_r_group_level",
            "spearman_rho_group_level",
            "within_5y_fraction_group_level",
            "baseline_mae_group_level",
            "delta_mae_vs_baseline_group_level",
            "n_samples",
            "n_groups",
            "n_donor_samples",
            "n_donors",
            "n_celltypes",
            "n_splits",
            "primary_level",
            "donor_balanced_training",
            "nested_cv_enabled",
            "inner_cv_folds",
            "selection_method",
            "bootstrap_unit",
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
    sampling_strategy = str(pcfg.get("sampling_strategy", "stratified")).lower()
    seed = int(cfg.get("run", {}).get("seed", 42))

    adata = ad.read_h5ad(inp_h5ad, backed="r")
    try:
        attach_biological_replicates(adata, cfg)
        attach_analysis_cell_types(adata, cfg)
        obs_cols = list(adata.obs.columns)

        age_col = pcfg.get("age_col")
        configured_sample_unit_col = pcfg.get("sample_unit_col")
        configured_group_col = pcfg.get("group_col")
        sample_unit_col = configured_sample_unit_col or pcfg.get("donor_col")
        group_col = configured_group_col or sample_unit_col
        celltype_col = pcfg.get("celltype_col")

        if age_col not in obs_cols:
            age_col = _first_existing(obs_cols, ["age"])
        if configured_sample_unit_col and sample_unit_col not in obs_cols:
            return (
                pd.DataFrame(),
                [],
                f"Configured sample-unit column is missing: {sample_unit_col}",
            )
        if sample_unit_col not in obs_cols:
            sample_unit_col = _first_existing(
                obs_cols, ["donor_id", "donor_id_y", "donor_id_x", "donor", "sample_id"]
            )
        if configured_group_col and group_col not in obs_cols:
            return (
                pd.DataFrame(),
                [],
                f"Configured prediction-group column is missing: {group_col}",
            )
        if group_col not in obs_cols:
            group_col = sample_unit_col
        if celltype_col not in obs_cols:
            celltype_col = _first_existing(
                obs_cols, ["cell_type", "majority_voting", "predicted_labels", "leiden"]
            )

        if age_col is None or sample_unit_col is None or group_col is None or celltype_col is None:
            return pd.DataFrame(), [], "Missing required obs columns for age prediction."

        if latent_key not in adata.obsm:
            return pd.DataFrame(), [], f"Missing latent representation in obsm: {latent_key}"

        obs = pd.DataFrame(
            {
                "age": adata.obs[age_col],
                "sample_unit_id": adata.obs[sample_unit_col],
                "donor_id": adata.obs[group_col],
                "cell_type": adata.obs[celltype_col],
            },
            index=adata.obs_names,
        )
        obs["age"] = pd.to_numeric(obs["age"], errors="coerce")

        valid_mask = (
            obs["age"].notna()
            & obs["sample_unit_id"].notna()
            & obs["donor_id"].notna()
            & obs["cell_type"].notna()
        )
        valid_idx = np.flatnonzero(valid_mask.to_numpy())
        if len(valid_idx) == 0:
            return pd.DataFrame(), [], "No valid cells after age/donor/cell_type filtering."

        valid_obs = obs.iloc[valid_idx]
        sample_contract = valid_obs.groupby("sample_unit_id", observed=False).agg(
            n_subjects=("donor_id", "nunique"),
            n_ages=("age", "nunique"),
        )
        invalid_samples = sample_contract[
            sample_contract["n_subjects"].ne(1) | sample_contract["n_ages"].ne(1)
        ]
        if not invalid_samples.empty:
            return (
                pd.DataFrame(),
                [],
                "Age prediction found sample units with conflicting subject or age assignments.",
            )

        if len(valid_idx) > max_cells:
            if sampling_strategy == "random":
                rng = np.random.default_rng(seed)
                valid_idx = np.sort(rng.choice(valid_idx, size=max_cells, replace=False))
            else:
                donor_values = obs.iloc[valid_idx]["sample_unit_id"].astype(str).to_numpy()
                celltype_values = obs.iloc[valid_idx]["cell_type"].astype(str).to_numpy()
                valid_idx = _stratified_sample_indices(
                    valid_idx=valid_idx,
                    donor_values=donor_values,
                    celltype_values=celltype_values,
                    max_cells=max_cells,
                    min_cells_per_group=min_cells_per_group,
                    seed=seed,
                )

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
            cell_df.groupby(["sample_unit_id", "donor_id", "cell_type"], observed=False)
            .agg(**agg_map)
            .reset_index()
        )
        agg = agg[agg["n_cells"] >= min_cells_per_group].copy()

        if agg.empty:
            return pd.DataFrame(), [], "No sample-unit-celltype groups pass min_cells_per_group."

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
        hcfg = pcfg.get(
            "hist_gbr_params", [{"max_depth": 4, "learning_rate": 0.05, "max_iter": 300}]
        )
        for i, p in enumerate(hcfg):
            params = dict(p)
            params.setdefault("random_state", seed)
            label = f"hist_gbr_{i + 1}"
            candidates.append(
                {
                    "model": "hist_gbr",
                    "family": "tree",
                    "label": label,
                    "params": params,
                    "build": lambda params=params: Pipeline(
                        steps=[
                            ("preproc", preproc_tree),
                            ("reg", HistGradientBoostingRegressor(**params)),
                        ]
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
                label = f"xgboost_{i + 1}"
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
        except ImportError as exc:
            raise ImportError(
                "age_prediction.model_order requests 'xgboost', but xgboost "
                "is not installed. Install the tracked environment or remove "
                "xgboost from model_order in an explicitly reviewed config."
            ) from exc

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
                label = f"lightgbm_{i + 1}"
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
        except ImportError as exc:
            raise ImportError(
                "age_prediction.model_order requests 'lightgbm', but lightgbm "
                "is not installed. Install it or remove lightgbm from "
                "model_order in an explicitly reviewed config."
            ) from exc

    return candidates


def _compute_donor_balanced_sample_weights(donor_ids: np.ndarray) -> np.ndarray:
    donors = pd.Series(np.asarray(donor_ids, dtype=str))
    counts = donors.value_counts()
    weights = donors.map(lambda d: 1.0 / float(counts.loc[d])).to_numpy(dtype=float)
    valid = np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return np.ones(donors.shape[0], dtype=float)
    mean_w = float(np.mean(weights[valid]))
    if not np.isfinite(mean_w) or mean_w <= 0.0:
        return np.ones(donors.shape[0], dtype=float)
    return weights / mean_w


def _fit_with_optional_sample_weight(
    model, x_train: pd.DataFrame, y_train: np.ndarray, sample_weight: np.ndarray | None
) -> None:
    if sample_weight is None:
        model.fit(x_train, y_train)
        return
    sw = np.asarray(sample_weight, dtype=float)
    if sw.shape[0] != y_train.shape[0]:
        model.fit(x_train, y_train)
        return

    for kwargs in ({"reg__sample_weight": sw}, {"sample_weight": sw}, {}):
        try:
            model.fit(x_train, y_train, **kwargs)
            return
        except (TypeError, ValueError):
            continue
    model.fit(x_train, y_train)


def _safe_correlation(corr_fn, x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(valid)) < 2:
        return np.nan
    x = x[valid]
    y = y[valid]
    if np.allclose(np.std(x), 0.0) or np.allclose(np.std(y), 0.0):
        return np.nan
    try:
        return float(corr_fn(x, y)[0])
    except Exception:
        return np.nan


def _compute_prediction_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_base: np.ndarray) -> dict:
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)
    yb = np.asarray(y_base, dtype=float)
    valid = np.isfinite(yt) & np.isfinite(yp) & np.isfinite(yb)
    if int(np.sum(valid)) == 0:
        return {
            "mae": np.nan,
            "median_ae": np.nan,
            "rmse": np.nan,
            "r2": np.nan,
            "pearson_r": np.nan,
            "spearman_rho": np.nan,
            "within_5y_fraction": np.nan,
            "baseline_mae": np.nan,
            "delta_mae_vs_baseline": np.nan,
            "n_samples": 0,
        }

    yt = yt[valid]
    yp = yp[valid]
    yb = yb[valid]

    abs_err = np.abs(yp - yt)
    base_abs_err = np.abs(yb - yt)
    ss_tot = float(np.sum((yt - float(np.mean(yt))) ** 2))
    if ss_tot <= 0.0:
        r2 = np.nan
    else:
        ss_res = float(np.sum((yt - yp) ** 2))
        r2 = 1.0 - (ss_res / ss_tot)

    mae = float(np.mean(abs_err))
    baseline_mae = float(np.mean(base_abs_err))
    return {
        "mae": mae,
        "median_ae": float(np.median(abs_err)),
        "rmse": float(np.sqrt(np.mean((yp - yt) ** 2))),
        "r2": float(r2),
        "pearson_r": _safe_correlation(stats.pearsonr, yt, yp),
        "spearman_rho": _safe_correlation(stats.spearmanr, yt, yp),
        "within_5y_fraction": float(np.mean(abs_err <= 5.0)),
        "baseline_mae": baseline_mae,
        "delta_mae_vs_baseline": float(baseline_mae - mae),
        "n_samples": int(yt.shape[0]),
    }


def _compute_r2_only(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)
    valid = np.isfinite(yt) & np.isfinite(yp)
    yt = yt[valid]
    yp = yp[valid]
    if yt.shape[0] < 2:
        return np.nan
    ss_tot = float(np.sum((yt - float(np.mean(yt))) ** 2))
    if ss_tot <= 0.0:
        return np.nan
    ss_res = float(np.sum((yt - yp) ** 2))
    return float(1.0 - (ss_res / ss_tot))


def _bootstrap_metric_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metric: str,
    n_boot: int,
    ci: float,
    seed: int,
    groups: np.ndarray | None = None,
) -> tuple[float, float]:
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)
    valid = np.isfinite(yt) & np.isfinite(yp)
    yt = yt[valid]
    yp = yp[valid]
    group_values = None
    if groups is not None:
        group_values = np.asarray(groups, dtype=str)
        if group_values.shape[0] != valid.shape[0]:
            raise ValueError("Bootstrap groups must match the prediction rows.")
        group_values = group_values[valid]
    n = int(yt.shape[0])
    if n == 0:
        return np.nan, np.nan

    def _eval(a: np.ndarray, b: np.ndarray) -> float:
        if metric == "mae":
            return float(np.mean(np.abs(a - b)))
        if metric == "r2":
            return _compute_r2_only(a, b)
        raise ValueError(f"Unsupported bootstrap metric: {metric}")

    if n == 1:
        v = _eval(yt, yp)
        return float(v), float(v)

    rng = np.random.default_rng(seed)
    n_iterations = int(max(n_boot, 1))
    boot_vals = np.empty(n_iterations, dtype=float)
    if group_values is None:
        indices = rng.integers(0, n, size=(n_iterations, n))
        for i in range(n_iterations):
            boot_vals[i] = _eval(yt[indices[i]], yp[indices[i]])
    else:
        unique_groups = np.unique(group_values)
        group_indices = {group: np.flatnonzero(group_values == group) for group in unique_groups}
        for i in range(n_iterations):
            sampled_groups = rng.choice(
                unique_groups,
                size=unique_groups.size,
                replace=True,
            )
            sampled_indices = np.concatenate([group_indices[group] for group in sampled_groups])
            boot_vals[i] = _eval(yt[sampled_indices], yp[sampled_indices])
    boot_vals = boot_vals[np.isfinite(boot_vals)]
    if boot_vals.size == 0:
        return np.nan, np.nan
    alpha = (1.0 - float(ci)) / 2.0
    low = float(np.quantile(boot_vals, alpha))
    high = float(np.quantile(boot_vals, 1.0 - alpha))
    return low, high


def _bootstrap_mae_r2_ci_from_pred_frame(
    pred_df: pd.DataFrame,
    n_boot: int,
    ci: float,
    seed: int,
) -> dict:
    if pred_df.empty:
        return {
            "mae_ci_low": np.nan,
            "mae_ci_high": np.nan,
            "r2_ci_low": np.nan,
            "r2_ci_high": np.nan,
        }
    y_true = pred_df["age_true"].to_numpy(dtype=float)
    y_pred = pred_df["age_pred"].to_numpy(dtype=float)
    groups = pred_df["donor_id"].astype(str).to_numpy() if "donor_id" in pred_df.columns else None
    mae_low, mae_high = _bootstrap_metric_ci(
        y_true,
        y_pred,
        metric="mae",
        n_boot=n_boot,
        ci=ci,
        seed=seed,
        groups=groups,
    )
    r2_low, r2_high = _bootstrap_metric_ci(
        y_true,
        y_pred,
        metric="r2",
        n_boot=n_boot,
        ci=ci,
        seed=seed + 17,
        groups=groups,
    )
    return {
        "mae_ci_low": mae_low,
        "mae_ci_high": mae_high,
        "r2_ci_low": r2_low,
        "r2_ci_high": r2_high,
    }


def _build_prediction_frames(
    agg_subset: pd.DataFrame,
    y_pred: np.ndarray,
    y_base: np.ndarray,
    fold_id: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = ["donor_id", "cell_type", "age", "n_cells"]
    if "sample_unit_id" in agg_subset.columns:
        columns.insert(1, "sample_unit_id")
    pred_group = agg_subset[columns].copy()
    pred_group = pred_group.rename(columns={"age": "age_true"})
    pred_group["evaluation_level"] = (
        "sample_unit_celltype" if "sample_unit_id" in pred_group.columns else "donor_celltype"
    )
    pred_group["age_pred"] = y_pred
    pred_group["age_base"] = y_base
    pred_group["residual"] = pred_group["age_pred"] - pred_group["age_true"]
    pred_group["abs_error"] = np.abs(pred_group["residual"])
    pred_group["within_5y"] = pred_group["abs_error"] <= 5.0
    pred_group["fold"] = fold_id
    pred_group["n_groups"] = 1

    pred_donor = _aggregate_group_predictions_to_donor(pred_group)
    return pred_group, pred_donor


def _run_group_cv_predictions_for_candidate(
    agg_subset: pd.DataFrame,
    latent_cols: list[str],
    splits: list[tuple[np.ndarray, np.ndarray]],
    cand: dict,
    clip_predictions: bool,
    donor_balanced_training: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = agg_subset["age"].to_numpy(dtype=float)
    groups = agg_subset["donor_id"].astype(str).to_numpy()
    x = agg_subset[latent_cols + ["cell_type"]].copy()

    y_pred = np.full_like(y, np.nan, dtype=float)
    y_base = np.full_like(y, np.nan, dtype=float)
    fold_id = np.full(y.shape[0], -1, dtype=int)

    for fold, (tr, te) in enumerate(splits):
        model = cand["build"]()
        train_weights = None
        if donor_balanced_training:
            train_weights = _compute_donor_balanced_sample_weights(groups[tr])
        _fit_with_optional_sample_weight(model, x.iloc[tr], y[tr], sample_weight=train_weights)
        fold_pred = model.predict(x.iloc[te])
        if clip_predictions:
            lo = float(np.min(y[tr]))
            hi = float(np.max(y[tr]))
            fold_pred = np.clip(fold_pred, lo, hi)
        y_pred[te] = fold_pred
        donor_train_age = agg_subset.iloc[tr].groupby("donor_id", observed=False)["age"].median()
        if donor_train_age.empty:
            y_base[te] = np.mean(y[tr])
        else:
            y_base[te] = float(donor_train_age.mean())
        fold_id[te] = fold

    return y_pred, y_base, fold_id


def _metrics_row_from_predictions(
    pred_group: pd.DataFrame,
    pred_donor: pd.DataFrame,
    *,
    model: str,
    family: str,
    params_json: str,
    alpha: float,
    primary_level: str,
    donor_balanced_training: bool,
    nested_cv_enabled: bool,
    inner_cv_folds: int,
    selection_method: str,
    n_donors: int,
    n_celltypes: int,
    n_splits: int,
    bootstrap_iterations: int,
    bootstrap_ci: float,
    seed: int,
) -> dict:
    group_metrics = _compute_prediction_metrics(
        y_true=pred_group["age_true"].to_numpy(dtype=float),
        y_pred=pred_group["age_pred"].to_numpy(dtype=float),
        y_base=pred_group["age_base"].to_numpy(dtype=float),
    )
    donor_metrics = _compute_prediction_metrics(
        y_true=pred_donor["age_true"].to_numpy(dtype=float),
        y_pred=pred_donor["age_pred"].to_numpy(dtype=float),
        y_base=pred_donor["age_base"].to_numpy(dtype=float),
    )
    if primary_level in {"donor", "sample_unit"}:
        primary_metrics = donor_metrics
        primary_pred = pred_donor
    else:
        primary_metrics = group_metrics
        primary_pred = pred_group

    primary_ci = _bootstrap_mae_r2_ci_from_pred_frame(
        primary_pred,
        n_boot=bootstrap_iterations,
        ci=bootstrap_ci,
        seed=seed,
    )
    donor_ci = (
        primary_ci
        if primary_level in {"donor", "sample_unit"}
        else _bootstrap_mae_r2_ci_from_pred_frame(
            pred_donor,
            n_boot=bootstrap_iterations,
            ci=bootstrap_ci,
            seed=seed + 97,
        )
    )

    return {
        "model": model,
        "family": family,
        "params": params_json,
        "alpha": alpha,
        "mae": primary_metrics["mae"],
        "mae_ci_low": primary_ci["mae_ci_low"],
        "mae_ci_high": primary_ci["mae_ci_high"],
        "median_ae": primary_metrics["median_ae"],
        "rmse": primary_metrics["rmse"],
        "r2": primary_metrics["r2"],
        "r2_ci_low": primary_ci["r2_ci_low"],
        "r2_ci_high": primary_ci["r2_ci_high"],
        "pearson_r": primary_metrics["pearson_r"],
        "spearman_rho": primary_metrics["spearman_rho"],
        "within_5y_fraction": primary_metrics["within_5y_fraction"],
        "baseline_mae": primary_metrics["baseline_mae"],
        "delta_mae_vs_baseline": primary_metrics["delta_mae_vs_baseline"],
        "mae_donor_level": donor_metrics["mae"],
        "mae_donor_ci_low": donor_ci["mae_ci_low"],
        "mae_donor_ci_high": donor_ci["mae_ci_high"],
        "median_ae_donor_level": donor_metrics["median_ae"],
        "rmse_donor_level": donor_metrics["rmse"],
        "r2_donor_level": donor_metrics["r2"],
        "r2_donor_ci_low": donor_ci["r2_ci_low"],
        "r2_donor_ci_high": donor_ci["r2_ci_high"],
        "pearson_r_donor_level": donor_metrics["pearson_r"],
        "spearman_rho_donor_level": donor_metrics["spearman_rho"],
        "within_5y_fraction_donor_level": donor_metrics["within_5y_fraction"],
        "baseline_mae_donor_level": donor_metrics["baseline_mae"],
        "delta_mae_vs_baseline_donor_level": donor_metrics["delta_mae_vs_baseline"],
        "mae_group_level": group_metrics["mae"],
        "median_ae_group_level": group_metrics["median_ae"],
        "rmse_group_level": group_metrics["rmse"],
        "r2_group_level": group_metrics["r2"],
        "pearson_r_group_level": group_metrics["pearson_r"],
        "spearman_rho_group_level": group_metrics["spearman_rho"],
        "within_5y_fraction_group_level": group_metrics["within_5y_fraction"],
        "baseline_mae_group_level": group_metrics["baseline_mae"],
        "delta_mae_vs_baseline_group_level": group_metrics["delta_mae_vs_baseline"],
        "n_samples": int(primary_metrics["n_samples"]),
        "n_groups": int(pred_group.shape[0]),
        "n_donor_samples": int(pred_donor.shape[0]),
        "n_donors": int(n_donors),
        "n_subjects": int(n_donors),
        "n_celltypes": int(n_celltypes),
        "n_splits": int(n_splits),
        "primary_level": primary_level,
        "donor_balanced_training": donor_balanced_training,
        "nested_cv_enabled": nested_cv_enabled,
        "inner_cv_folds": int(inner_cv_folds),
        "selection_method": selection_method,
        "bootstrap_unit": "subject",
    }


def _make_group_splits(
    groups: np.ndarray,
    *,
    requested_folds: int,
    context: str,
) -> list[tuple[np.ndarray, np.ndarray]]:
    from sklearn.model_selection import GroupKFold

    groups = np.asarray(groups, dtype=str)
    unique_groups = np.unique(groups)
    n_splits = min(int(requested_folds), int(unique_groups.size))
    if n_splits < 2:
        raise ValueError(f"{context} requires at least two distinct subject groups.")

    splitter = GroupKFold(n_splits=n_splits)
    splits = list(splitter.split(np.zeros(groups.shape[0]), groups=groups))
    test_fold = np.full(groups.shape[0], -1, dtype=int)
    for fold, (train_idx, test_idx) in enumerate(splits):
        overlap = set(groups[train_idx]).intersection(groups[test_idx])
        if overlap:
            raise ValueError(f"{context} fold {fold} leaks {len(overlap)} subject(s).")
        test_fold[test_idx] = fold

    if bool(np.any(test_fold < 0)):
        raise RuntimeError(f"{context} did not assign every row to one test fold.")
    fold_counts = (
        pd.DataFrame({"subject_id": groups, "fold": test_fold})
        .groupby("subject_id", observed=False)["fold"]
        .nunique()
    )
    if bool(fold_counts.ne(1).any()):
        raise RuntimeError(f"{context} assigned a subject to more than one test fold.")
    return splits


def _select_inner_best_candidate(
    agg_train: pd.DataFrame,
    latent_cols: list[str],
    candidates: list[dict],
    inner_cv_folds: int,
    clip_predictions: bool,
    donor_balanced_training: bool,
    primary_level: str,
) -> tuple[dict, float]:
    groups_train = agg_train["donor_id"].astype(str).to_numpy()
    unique_train_groups = np.unique(groups_train)
    if unique_train_groups.size < 2:
        return candidates[0], np.inf
    inner_splits = _make_group_splits(
        groups_train,
        requested_folds=inner_cv_folds,
        context="Inner cross-validation",
    )

    best_cand = candidates[0]
    best_score = np.inf
    for cand in candidates:
        y_pred, y_base, fold_id = _run_group_cv_predictions_for_candidate(
            agg_subset=agg_train,
            latent_cols=latent_cols,
            splits=inner_splits,
            cand=cand,
            clip_predictions=clip_predictions,
            donor_balanced_training=donor_balanced_training,
        )
        pred_group, pred_donor = _build_prediction_frames(
            agg_subset=agg_train, y_pred=y_pred, y_base=y_base, fold_id=fold_id
        )
        score_metrics = (
            _compute_prediction_metrics(
                y_true=pred_donor["age_true"].to_numpy(dtype=float),
                y_pred=pred_donor["age_pred"].to_numpy(dtype=float),
                y_base=pred_donor["age_base"].to_numpy(dtype=float),
            )
            if primary_level in {"donor", "sample_unit"}
            else _compute_prediction_metrics(
                y_true=pred_group["age_true"].to_numpy(dtype=float),
                y_pred=pred_group["age_pred"].to_numpy(dtype=float),
                y_base=pred_group["age_base"].to_numpy(dtype=float),
            )
        )
        score = float(score_metrics["mae"])
        if np.isfinite(score) and score < best_score:
            best_score = score
            best_cand = cand
    return best_cand, best_score


def _aggregate_group_predictions_to_donor(pred_group: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    aggregate_col = "sample_unit_id" if "sample_unit_id" in pred_group.columns else "donor_id"
    for aggregate_id, sub in pred_group.groupby(aggregate_col, observed=True):
        if sub.empty:
            continue
        weights = sub["n_cells"].to_numpy(dtype=float)
        weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
        if float(weights.sum()) <= 0.0:
            weights = np.ones(sub.shape[0], dtype=float)
        age_true = float(np.median(sub["age_true"].to_numpy(dtype=float)))
        age_pred = float(np.average(sub["age_pred"].to_numpy(dtype=float), weights=weights))
        age_base = float(np.average(sub["age_base"].to_numpy(dtype=float), weights=weights))
        rows.append(
            {
                "evaluation_level": (
                    "sample_unit" if aggregate_col == "sample_unit_id" else "donor"
                ),
                "donor_id": str(sub["donor_id"].iloc[0]),
                **(
                    {"sample_unit_id": str(aggregate_id)}
                    if aggregate_col == "sample_unit_id"
                    else {}
                ),
                "cell_type": "__all__",
                "age_true": age_true,
                "age_pred": age_pred,
                "age_base": age_base,
                "n_cells": int(np.sum(sub["n_cells"].to_numpy(dtype=int))),
                "n_groups": int(sub.shape[0]),
                "fold": int(sub["fold"].iloc[0]) if sub.shape[0] > 0 else -1,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["residual"] = out["age_pred"] - out["age_true"]
    out["abs_error"] = np.abs(out["residual"])
    out["within_5y"] = out["abs_error"] <= 5.0
    return out


def _evaluate_candidates(
    agg: pd.DataFrame, latent_cols: list[str], cfg: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pcfg = cfg.get("age_prediction", {})
    cv_folds = int(pcfg.get("cv_folds", 5))
    nested_cv_enabled = bool(pcfg.get("nested_cv", True))
    inner_cv_folds = int(pcfg.get("nested_inner_cv_folds", 3))
    clip_predictions = bool(pcfg.get("clip_predictions_to_train_age_range", True))
    donor_balanced_training = bool(pcfg.get("donor_balanced_training", True))
    default_level = "sample_unit" if "sample_unit_id" in agg.columns else "donor"
    primary_level_cfg = str(pcfg.get("primary_metric_level", default_level)).strip().lower()
    allowed_levels = {"donor", "donor_celltype", "sample_unit", "sample_unit_celltype"}
    primary_level = default_level if primary_level_cfg not in allowed_levels else primary_level_cfg
    bootstrap_iterations = int(pcfg.get("bootstrap_iterations", 2000))
    bootstrap_ci = float(pcfg.get("bootstrap_ci", 0.95))
    seed = int(cfg.get("run", {}).get("seed", 42))

    y = agg["age"].to_numpy(dtype=float)
    groups = agg["donor_id"].astype(str).to_numpy()
    x = agg[latent_cols + ["cell_type"]].copy()

    unique_groups = np.unique(groups)
    if unique_groups.size < 3:
        raise ValueError("Need at least 3 donors for cross-donor CV.")
    splits = _make_group_splits(
        groups,
        requested_folds=cv_folds,
        context="Outer cross-validation",
    )
    n_splits = len(splits)

    candidates = _load_model_builders(cfg, latent_cols=latent_cols, seed=seed)
    if not candidates:
        raise ValueError("No valid models available (check model_order and installed libraries).")

    metrics_rows = []
    pred_rows = []
    for i, cand in enumerate(candidates):
        y_pred, y_base, fold_id = _run_group_cv_predictions_for_candidate(
            agg_subset=agg,
            latent_cols=latent_cols,
            splits=splits,
            cand=cand,
            clip_predictions=clip_predictions,
            donor_balanced_training=donor_balanced_training,
        )
        pred_group, pred_donor = _build_prediction_frames(
            agg_subset=agg, y_pred=y_pred, y_base=y_base, fold_id=fold_id
        )
        params_json = json.dumps(cand["params"], sort_keys=True)
        metrics_rows.append(
            _metrics_row_from_predictions(
                pred_group=pred_group,
                pred_donor=pred_donor,
                model=cand["model"],
                family=cand["family"],
                params_json=params_json,
                alpha=float(cand["params"].get("alpha", np.nan)),
                primary_level=primary_level,
                donor_balanced_training=donor_balanced_training,
                nested_cv_enabled=False,
                inner_cv_folds=0,
                selection_method="fixed_outer_cv",
                n_donors=int(unique_groups.size),
                n_celltypes=int(agg["cell_type"].nunique()),
                n_splits=int(n_splits),
                bootstrap_iterations=bootstrap_iterations,
                bootstrap_ci=bootstrap_ci,
                seed=seed + (i * 131),
            )
        )

        pred_group["model"] = cand["model"]
        pred_group["params"] = params_json
        pred_group["selected_model"] = np.nan
        pred_group["selected_params"] = np.nan
        pred_donor["model"] = cand["model"]
        pred_donor["params"] = params_json
        pred_donor["selected_model"] = np.nan
        pred_donor["selected_params"] = np.nan
        pred_rows.append(pred_group)
        pred_rows.append(pred_donor)

    if nested_cv_enabled:
        y_pred_nested = np.full_like(y, np.nan, dtype=float)
        y_base_nested = np.full_like(y, np.nan, dtype=float)
        fold_nested = np.full(y.shape[0], -1, dtype=int)
        selected_model = np.empty(y.shape[0], dtype=object)
        selected_params = np.empty(y.shape[0], dtype=object)
        selected_model[:] = ""
        selected_params[:] = ""
        fold_selection_records: list[dict] = []

        for outer_fold, (tr, te) in enumerate(splits):
            agg_train = agg.iloc[tr].reset_index(drop=True)
            best_cand, inner_score = _select_inner_best_candidate(
                agg_train=agg_train,
                latent_cols=latent_cols,
                candidates=candidates,
                inner_cv_folds=inner_cv_folds,
                clip_predictions=clip_predictions,
                donor_balanced_training=donor_balanced_training,
                primary_level=primary_level,
            )
            model = best_cand["build"]()
            train_weights = None
            if donor_balanced_training:
                train_weights = _compute_donor_balanced_sample_weights(groups[tr])
            _fit_with_optional_sample_weight(model, x.iloc[tr], y[tr], sample_weight=train_weights)
            fold_pred = model.predict(x.iloc[te])
            if clip_predictions:
                lo = float(np.min(y[tr]))
                hi = float(np.max(y[tr]))
                fold_pred = np.clip(fold_pred, lo, hi)
            y_pred_nested[te] = fold_pred
            donor_train_age = agg.iloc[tr].groupby("donor_id", observed=False)["age"].median()
            if donor_train_age.empty:
                y_base_nested[te] = np.mean(y[tr])
            else:
                y_base_nested[te] = float(donor_train_age.mean())
            fold_nested[te] = outer_fold

            best_params_json = json.dumps(best_cand["params"], sort_keys=True)
            selected_model[te] = str(best_cand["model"])
            selected_params[te] = best_params_json
            fold_selection_records.append(
                {
                    "outer_fold": int(outer_fold),
                    "selected_model": str(best_cand["model"]),
                    "selected_params": best_params_json,
                    "inner_mae": float(inner_score),
                }
            )

        pred_group_nested, pred_donor_nested = _build_prediction_frames(
            agg_subset=agg,
            y_pred=y_pred_nested,
            y_base=y_base_nested,
            fold_id=fold_nested,
        )
        pred_group_nested["selected_model"] = selected_model
        pred_group_nested["selected_params"] = selected_params
        fold_to_model = {r["outer_fold"]: r["selected_model"] for r in fold_selection_records}
        fold_to_params = {r["outer_fold"]: r["selected_params"] for r in fold_selection_records}
        pred_donor_nested["selected_model"] = (
            pred_donor_nested["fold"].map(fold_to_model).fillna("")
        )
        pred_donor_nested["selected_params"] = (
            pred_donor_nested["fold"].map(fold_to_params).fillna("")
        )

        nested_params_json = json.dumps(
            {
                "inner_cv_folds": int(inner_cv_folds),
                "selection_metric": "mae",
                "selection_level": primary_level,
            },
            sort_keys=True,
        )
        metrics_rows.append(
            _metrics_row_from_predictions(
                pred_group=pred_group_nested,
                pred_donor=pred_donor_nested,
                model="nested_selected",
                family="meta",
                params_json=nested_params_json,
                alpha=np.nan,
                primary_level=primary_level,
                donor_balanced_training=donor_balanced_training,
                nested_cv_enabled=True,
                inner_cv_folds=int(inner_cv_folds),
                selection_method="nested_outer_inner",
                n_donors=int(unique_groups.size),
                n_celltypes=int(agg["cell_type"].nunique()),
                n_splits=int(n_splits),
                bootstrap_iterations=bootstrap_iterations,
                bootstrap_ci=bootstrap_ci,
                seed=seed + 7001,
            )
        )
        pred_group_nested["model"] = "nested_selected"
        pred_group_nested["params"] = nested_params_json
        pred_donor_nested["model"] = "nested_selected"
        pred_donor_nested["params"] = nested_params_json
        pred_rows.append(pred_group_nested)
        pred_rows.append(pred_donor_nested)

    metrics = pd.DataFrame(metrics_rows).sort_values("mae", ascending=True).reset_index(drop=True)
    metrics["is_best"] = False
    if not metrics.empty:
        nested_mask = metrics["model"].eq("nested_selected")
        if nested_cv_enabled and nested_mask.any():
            metrics.loc[nested_mask, "is_best"] = True
        else:
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
    best_mask = (pred_all["model"] == best_row["model"]) & (
        pred_all["params"] == best_row["params"]
    )
    pred_all["is_best_model"] = best_mask
    pred_all = pred_all.sort_values(
        ["is_best_model", "model", "evaluation_level", "donor_id", "cell_type"],
        ascending=[False, True, True, True, True],
    )

    return pred_all, metrics


def _bootstrap_mean_ci(
    values: np.ndarray, n_boot: int, ci: float, seed: int
) -> tuple[float, float, float]:
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


def _build_model_comparison_summary(
    pred_all: pd.DataFrame, metrics: pd.DataFrame, cfg: dict
) -> pd.DataFrame:
    pcfg = cfg.get("age_prediction", {})
    n_boot = int(pcfg.get("bootstrap_iterations", 2000))
    ci = float(pcfg.get("bootstrap_ci", 0.95))
    seed = int(cfg.get("run", {}).get("seed", 42))

    if pred_all.empty or metrics.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    pred_all = pred_all.copy()
    if "evaluation_level" in pred_all.columns:
        primary_rows = pred_all[pred_all["evaluation_level"].astype(str) == "donor"].copy()
        if primary_rows.empty:
            primary_rows = pred_all[
                pred_all["evaluation_level"].astype(str) == "donor_celltype"
            ].copy()
        if primary_rows.empty:
            primary_rows = pred_all.copy()
    else:
        primary_rows = pred_all.copy()
    primary_rows["model_key"] = primary_rows["model"] + "||" + primary_rows["params"].astype(str)

    fold_mae = (
        primary_rows.groupby(["model", "params", "model_key", "fold"], observed=False)["abs_error"]
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
        ci_low, ci_high, prob_gt0 = _bootstrap_mean_ci(
            delta, n_boot=n_boot, ci=ci, seed=seed + len(summary_rows)
        )

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
    return out.sort_values(
        ["candidate_is_best", "delta_mae_mean"], ascending=[False, False]
    ).reset_index(drop=True)


def _plot_observed_vs_pred(
    pred_best: pd.DataFrame, metrics: pd.DataFrame, path: Path, dpi: int, density_label: str
) -> None:
    if pred_best.empty:
        _save_placeholder(
            path, "Age Prediction: Observed vs Predicted", "No prediction rows available.", dpi
        )
        return

    best = metrics.loc[metrics["is_best"]].iloc[0]
    x = pred_best["age_true"].to_numpy(dtype=float)
    y = pred_best["age_pred"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    hb = ax.hexbin(x, y, gridsize=32, cmap="Blues", mincnt=1)
    cbar = fig.colorbar(hb, ax=ax, shrink=0.8)
    cbar.set_label(density_label)
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
    mae_txt = f"{best['mae']:.2f}"
    if np.isfinite(float(best.get("mae_ci_low", np.nan))) and np.isfinite(
        float(best.get("mae_ci_high", np.nan))
    ):
        mae_txt += f" [{float(best['mae_ci_low']):.2f}, {float(best['mae_ci_high']):.2f}]"
    r2_txt = f"{best['r2']:.3f}"
    if np.isfinite(float(best.get("r2_ci_low", np.nan))) and np.isfinite(
        float(best.get("r2_ci_high", np.nan))
    ):
        r2_txt += f" [{float(best['r2_ci_low']):.3f}, {float(best['r2_ci_high']):.3f}]"
    txt = (
        f"Model={best['model']}\n"
        f"MAE={mae_txt}, median AE={best['median_ae']:.2f}, RMSE={best['rmse']:.2f}\n"
        f"R2={r2_txt}, Pearson r={best['pearson_r']:.3f}, within 5y={best['within_5y_fraction']:.2%}\n"
        f"Top MAE: {top_txt}"
    )
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=10)
    style_axis(ax, grid="y")
    ax.legend(frameon=False, loc="lower right")
    finalize_and_save(fig, path, dpi)


def _plot_mae_by_celltype(
    pred_best: pd.DataFrame,
    best_row: pd.Series,
    path: Path,
    top_n: int,
    min_groups_for_plot: int,
    dpi: int,
) -> None:
    ct_all = (
        pred_best.groupby("cell_type", observed=False)
        .agg(
            mae=("abs_error", "mean"),
            median_ae=("abs_error", "median"),
            n_groups=("abs_error", "size"),
            within_5y_fraction=("within_5y", "mean"),
        )
        .sort_values("mae", ascending=False)
    )
    if ct_all.empty:
        _save_placeholder(
            path, "Age Prediction: MAE by Cell Type", "No prediction rows available.", dpi
        )
        return

    ct = ct_all[ct_all["n_groups"] >= int(min_groups_for_plot)].copy()
    if ct.empty:
        _save_placeholder(
            path,
            "Age Prediction: MAE by Cell Type",
            f"No cell types meet n_groups >= {int(min_groups_for_plot)}.",
            dpi,
        )
        return

    ct = ct.head(top_n).sort_values("mae", ascending=True)
    fig, ax = plt.subplots(figsize=(10.5, 6.3))
    ax.barh(ct.index, ct["mae"].to_numpy(), color=PALETTE["primary"])
    ax.set_xlabel("Mean absolute error (years)")
    ax.set_ylabel("Cell type")
    ax.set_title(
        f"Cross-validated age prediction error by cell type\n"
        f"Best model: {best_row['model']} (n_groups>={int(min_groups_for_plot)})"
    )
    xmax = float(np.nanmax(ct["mae"].to_numpy())) if len(ct) > 0 else 1.0
    annotation_x = xmax * 1.03
    ax.set_xlim(0.0, xmax * 1.70)
    for i, (_, row) in enumerate(ct.iterrows()):
        ax.text(
            annotation_x,
            i,
            f"MAE={row['mae']:.1f}, n={int(row['n_groups'])}, medAE={row['median_ae']:.1f}, within5y={row['within_5y_fraction']:.1%}",
            va="center",
            ha="left",
            fontsize=8,
            color="#333333",
            clip_on=False,
        )
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
    min_groups_for_mae_plot = int(pcfg.get("min_groups_for_mae_plot", 20))

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
        require_placeholder_permission(cfg, "Age prediction requires scikit-learn.")
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
        require_placeholder_permission(
            cfg, err or "Age prediction has no aggregated replicate-cell-type rows."
        )
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
    print(
        "[age_prediction] "
        f"prepared {len(agg):,} sample-unit-cell-type groups from "
        f"{agg['donor_id'].nunique():,} subject groups",
        flush=True,
    )

    try:
        print("[age_prediction] evaluating grouped CV models", flush=True)
        pred_all, metrics = _evaluate_candidates(agg=agg, latent_cols=latent_cols, cfg=cfg)
    except Exception as exc:
        print(
            f"[age_prediction] grouped CV failed: {type(exc).__name__}: {exc}",
            flush=True,
        )
        require_placeholder_permission(
            cfg, f"Age-prediction grouped cross-validation failed: {exc}"
        )
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
    grouping_col = str(pcfg.get("group_col", pcfg.get("donor_col", "donor_id")))
    sample_unit_col = str(pcfg.get("sample_unit_col", pcfg.get("donor_col", "donor_id")))
    if grouping_col == sample_unit_col:
        pred_out = pred_all.drop(columns=["sample_unit_id"], errors="ignore").rename(
            columns={"donor_id": grouping_col}
        )
    else:
        pred_out = pred_all.rename(
            columns={"donor_id": grouping_col, "sample_unit_id": sample_unit_col}
        )
    metrics["grouping_id_column"] = grouping_col
    metrics["sample_unit_id_column"] = sample_unit_col
    summary["grouping_id_column"] = grouping_col
    summary["sample_unit_id_column"] = sample_unit_col
    pred_out.to_csv(table_pred, index=False)
    metrics.to_csv(table_metrics, index=False)
    summary.to_csv(table_summary, index=False)
    best_row = metrics.loc[metrics["is_best"]].iloc[0]
    pred_best = pred_all[pred_all["is_best_model"]].copy()
    if "evaluation_level" in pred_best.columns:
        pred_best_donor = pred_best[
            pred_best["evaluation_level"].astype(str).isin(["sample_unit", "donor"])
        ].copy()
        pred_best_group = pred_best[
            pred_best["evaluation_level"]
            .astype(str)
            .isin(["sample_unit_celltype", "donor_celltype"])
        ].copy()
    else:
        pred_best_donor = pd.DataFrame()
        pred_best_group = pd.DataFrame()
    if pred_best_donor.empty:
        pred_best_donor = pred_best
        density_label = "Sample-unit-cell-type groups"
    else:
        density_label = "Sample units"
    if pred_best_group.empty:
        pred_best_group = pred_best
    _plot_observed_vs_pred(pred_best_donor, metrics, fig_pred, dpi, density_label=density_label)
    _plot_mae_by_celltype(
        pred_best_group,
        best_row=best_row,
        path=fig_mae,
        top_n=top_n_celltypes,
        min_groups_for_plot=min_groups_for_mae_plot,
        dpi=dpi,
    )


if __name__ == "__main__":
    main()
