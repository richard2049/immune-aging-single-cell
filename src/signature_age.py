from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .biological_replicates import attach_biological_replicates
from .plot_style import (
    DIVERGING_CMAP,
    PALETTE,
    apply_publication_style,
    finalize_and_save,
    save_placeholder,
    style_axis,
)
from .scientific_guardrails import (
    build_complete_nuisance_design,
    require_placeholder_permission,
    residualize_complete,
    resolve_covariates,
    validate_replicate_covariates,
)
from .utils import ensure_dir, load_config

DEFAULT_SIGNATURES = {
    "ifn_response": [
        "ISG15",
        "IFIT1",
        "IFIT3",
        "IFI6",
        "IFI44L",
        "MX1",
        "OAS1",
        "STAT1",
        "IRF7",
        "CXCL10",
    ],
    "inflammatory_nfkb": [
        "NFKBIA",
        "TNF",
        "IL1B",
        "CXCL8",
        "PTGS2",
        "CCL2",
        "CCL3",
        "CCL4",
        "JUNB",
        "FOS",
    ],
    "mitochondrial_stress": [
        "MT-ND1",
        "MT-ND2",
        "MT-CO1",
        "MT-CO2",
        "MT-CYB",
        "HSPD1",
        "HSPE1",
        "LONP1",
    ],
    "proteostasis_upr": [
        "HSPA1A",
        "HSPA1B",
        "HSP90AA1",
        "DNAJB1",
        "XBP1",
        "ATF4",
        "DDIT3",
        "HSPH1",
        "BAG3",
    ],
    "sasp_proxy": [
        "IL6",
        "IL1A",
        "IL1B",
        "CXCL8",
        "CCL2",
        "MMP9",
        "TNF",
        "SERPINE1",
        "IGFBP7",
    ],
}

SIGNATURE_LABEL_OVERRIDES = {
    "ifn_response": "IFN response",
    "inflammatory_nfkb": "Inflammatory (NF-kB)",
    "mitochondrial_stress": "Mitochondrial stress",
    "proteostasis_upr": "Proteostasis (UPR)",
    "sasp_proxy": "SASP proxy",
}


def _save_placeholder(path: Path, title: str, message: str, dpi: int) -> None:
    save_placeholder(path=path, title=title, message=message, dpi=dpi)


def _format_signature_label(name: str) -> str:
    label = SIGNATURE_LABEL_OVERRIDES.get(name)
    if label:
        return label
    return name.replace("_", " ").strip().title()


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


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    m = pvals.size
    out = np.full(m, np.nan, dtype=float)
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


def _extract_signatures(cfg: dict) -> dict[str, list[str]]:
    scfg = cfg.get("signature_age", {})
    raw = scfg.get("signatures")
    if not isinstance(raw, dict) or not raw:
        return DEFAULT_SIGNATURES

    out: dict[str, list[str]] = {}
    for name, genes in raw.items():
        if not isinstance(name, str):
            continue
        if isinstance(genes, list):
            clean = [str(g).strip() for g in genes if str(g).strip()]
            if clean:
                out[name] = clean
    return out if out else DEFAULT_SIGNATURES


def _pick_columns(obs_cols: list[str], cfg: dict) -> dict[str, str | None]:
    scfg = cfg.get("signature_age", {})
    age = scfg.get("age_col")
    donor = scfg.get("donor_col")
    celltype = scfg.get("celltype_col")
    sex = scfg.get("sex_col")

    if age not in obs_cols:
        age = _first_existing(obs_cols, ["age"])
    if donor not in obs_cols:
        donor = _first_existing(
            obs_cols, ["donor_id", "donor_id_y", "donor_id_x", "donor", "sample_id"]
        )
    if celltype not in obs_cols:
        celltype = _first_existing(
            obs_cols, ["cell_type", "majority_voting", "predicted_labels", "leiden"]
        )
    if sex not in obs_cols:
        sex = _first_existing(obs_cols, ["sex"])

    return {"age": age, "donor": donor, "cell_type": celltype, "sex": sex}


def _subset_cells(
    adata,
    age_col: str,
    donor_col: str,
    celltype_col: str,
    max_cells: int,
    min_cells_per_group: int,
    seed: int,
    sampling_strategy: str = "stratified",
) -> np.ndarray:
    obs = adata.obs
    keep = (
        pd.to_numeric(obs[age_col], errors="coerce").notna()
        & obs[donor_col].notna()
        & obs[celltype_col].notna()
    )
    idx = np.flatnonzero(keep.to_numpy())
    if len(idx) <= max_cells:
        return idx

    sampling_strategy = str(sampling_strategy).lower()
    if sampling_strategy == "random":
        rng = np.random.default_rng(seed)
        chosen = rng.choice(idx, size=max_cells, replace=False)
        chosen.sort()
        return chosen

    donor_values = obs.iloc[idx][donor_col].astype(str).to_numpy()
    celltype_values = obs.iloc[idx][celltype_col].astype(str).to_numpy()
    return _stratified_sample_indices(
        valid_idx=idx,
        donor_values=donor_values,
        celltype_values=celltype_values,
        max_cells=max_cells,
        min_cells_per_group=min_cells_per_group,
        seed=seed,
    )


def _materialize_signature_genes(
    adata_backed,
    selected_indices: np.ndarray,
    signatures: dict[str, list[str]],
    assume_log1p: bool,
    target_sum: float,
    chunk_size: int,
) -> ad.AnnData:
    present_genes = sorted(
        {gene for genes in signatures.values() for gene in genes if gene in adata_backed.var_names}
    )
    if not present_genes:
        return ad.AnnData(
            X=sp.csr_matrix((len(selected_indices), 0), dtype=np.float32),
            obs=adata_backed.obs.iloc[selected_indices].copy(),
            var=pd.DataFrame(index=pd.Index([], dtype=str)),
        )

    gene_indices = adata_backed.var_names.get_indexer(present_genes)
    matrix_parts: list[sp.csr_matrix] = []
    obs_parts: list[pd.DataFrame] = []
    chunk_size = max(int(chunk_size), 1)
    for start in range(0, len(selected_indices), chunk_size):
        end = min(start + chunk_size, len(selected_indices))
        row_indices = selected_indices[start:end]
        chunk = adata_backed[row_indices].to_memory()
        obs_parts.append(chunk.obs.copy())

        full_matrix = chunk.X
        signature_matrix = full_matrix[:, gene_indices]
        if sp.issparse(signature_matrix):
            signature_matrix = signature_matrix.tocsr().astype(
                np.float32,
                copy=True,
            )
        else:
            signature_matrix = sp.csr_matrix(np.asarray(signature_matrix, dtype=np.float32))

        if not assume_log1p:
            totals = np.asarray(full_matrix.sum(axis=1)).ravel()
            scales = np.divide(
                float(target_sum),
                totals,
                out=np.zeros_like(totals, dtype=np.float64),
                where=totals > 0,
            )
            signature_matrix = sp.diags(scales).dot(signature_matrix).tocsr()
            signature_matrix.data = np.log1p(signature_matrix.data)
        matrix_parts.append(signature_matrix)
        print(
            f"[signature_age] materialized {end:,}/{len(selected_indices):,} selected cells",
            flush=True,
        )

    matrix = sp.vstack(matrix_parts, format="csr")
    obs = pd.concat(obs_parts, axis=0)
    var = adata_backed.var.iloc[gene_indices].copy()
    return ad.AnnData(X=matrix, obs=obs, var=var)


def _mode_or_na(values: pd.Series) -> str:
    mode = values.mode(dropna=True)
    if mode.empty:
        return "NA"
    return str(mode.iloc[0])


def _resolve_covariate_cols(
    obs_cols: list[str], cfg: dict, colmap: dict[str, str | None]
) -> list[str]:
    scfg = cfg.get("signature_age", {})
    defaults = (str(colmap["sex"]),) if colmap.get("sex") is not None else ()
    return resolve_covariates(
        obs_cols,
        scfg,
        section_name="signature_age",
        defaults=defaults,
    )


def _prepare_xy_for_stats(
    g: pd.DataFrame,
    score_col: str,
    covariate_cols: list[str],
    adjust_covariates: bool,
) -> tuple[np.ndarray, np.ndarray, bool]:
    x = g["age"].to_numpy(dtype=float)
    y = g[score_col].to_numpy(dtype=float)
    if not adjust_covariates or not covariate_cols:
        return x, y, False

    context = f"signature association for {score_col!r}"
    design = build_complete_nuisance_design(g[covariate_cols].copy(), context=context)
    x_res = residualize_complete(x, design, label="age", context=context)
    y_res = residualize_complete(y, design, label=score_col, context=context)
    return x_res, y_res, True


def _deduplicate_donors(
    df: pd.DataFrame, score_col: str, covariate_cols: list[str]
) -> pd.DataFrame:
    cols = ["donor_id", "age", score_col] + [c for c in covariate_cols if c in df.columns]
    out = df[cols].dropna(subset=["donor_id", "age", score_col]).copy()
    if out.empty:
        return out
    out = out.sort_values(["donor_id", "age"]).drop_duplicates(subset=["donor_id"], keep="last")
    return out


def _filter_assoc_for_panels(
    assoc: pd.DataFrame,
    significant_only: bool,
    min_donors_for_panel: int,
) -> pd.DataFrame:
    out = assoc.copy()
    if significant_only:
        out = out[out["fdr_significant"]].copy()
    if min_donors_for_panel > 0:
        out = out[out["n_donors"] >= int(min_donors_for_panel)].copy()
    return out


def _sort_assoc_for_panels(assoc: pd.DataFrame) -> pd.DataFrame:
    return assoc.sort_values(["fdr", "pvalue", "spearman_rho"], ascending=[True, True, False])


def _mean_expression_per_signature(
    adata, signatures: dict[str, list[str]], min_genes_present: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    var_upper = pd.Index([str(v).upper() for v in adata.var_names])
    rows = []
    score_df = pd.DataFrame(index=adata.obs_names)

    for sig_name, genes in signatures.items():
        gene_upper = [g.upper() for g in genes]
        mask = var_upper.isin(gene_upper)
        present = adata.var_names[mask].tolist()
        missing = [g for g in genes if g.upper() not in set(var_upper[mask])]

        rows.append(
            {
                "signature": sig_name,
                "n_genes_requested": len(genes),
                "n_genes_present": len(present),
                "present_genes": ";".join(map(str, present)),
                "missing_genes": ";".join(map(str, missing)),
            }
        )

        if len(present) < min_genes_present:
            score_df[f"score__{sig_name}"] = np.nan
            continue

        x = adata[:, present].X
        if sp.issparse(x):
            vals = np.asarray(x.mean(axis=1)).ravel()
        else:
            vals = np.mean(np.asarray(x), axis=1)
        score_df[f"score__{sig_name}"] = vals

    sig_meta = pd.DataFrame(rows)
    return score_df, sig_meta


def _aggregate_scores(
    obs: pd.DataFrame,
    score_cols: list[str],
    min_cells_per_group: int,
    covariate_cols: list[str],
) -> pd.DataFrame:
    agg_map: dict[str, tuple[str, str | callable]] = {
        "age": ("age", "median"),
        "n_cells": ("cell_type", "size"),
    }
    for c in covariate_cols:
        if c in obs.columns:
            agg_map[c] = (c, _mode_or_na)
    for c in score_cols:
        agg_map[c] = (c, "mean")

    grouped = obs.groupby(["donor_id", "cell_type"], observed=False).agg(**agg_map).reset_index()
    grouped = grouped[grouped["n_cells"] >= min_cells_per_group].copy()
    return grouped


def _empty_association_table() -> pd.DataFrame:
    cols = [
        "cell_type",
        "signature",
        "n_donors",
        "age_min",
        "age_max",
        "age_span",
        "spearman_rho",
        "pvalue",
        "spearman_rho_ci_low",
        "spearman_rho_ci_high",
        "slope_per_year",
        "slope_per_year_ci_low",
        "slope_per_year_ci_high",
        "effect_per_10y",
        "effect_per_10y_ci_low",
        "effect_per_10y_ci_high",
        "slope_pvalue",
        "r_squared",
        "adjusted",
        "adjusted_for",
        "fdr",
        "direction",
        "fdr_significant",
        "rank",
    ]
    return pd.DataFrame(columns=cols)


def _association_table(
    agg: pd.DataFrame,
    score_cols: list[str],
    min_donors_per_celltype: int,
    min_age_span: float,
    covariate_cols: list[str],
    adjust_covariates: bool,
    bootstrap_iterations: int,
    bootstrap_ci: float,
    seed: int,
) -> pd.DataFrame:
    rows = []
    adjusted_for = ",".join(covariate_cols) if (adjust_covariates and covariate_cols) else ""
    pair_idx = 0
    for cell_type, gct in agg.groupby("cell_type", observed=False):
        for sc_col in score_cols:
            t = _deduplicate_donors(gct, sc_col, covariate_cols=covariate_cols)
            n_donors = int(t["donor_id"].nunique())
            if n_donors < min_donors_per_celltype:
                continue

            age_span = float(t["age"].max() - t["age"].min())
            if age_span < min_age_span:
                continue

            x, y, used_adjustment = _prepare_xy_for_stats(
                t,
                score_col=sc_col,
                covariate_cols=covariate_cols,
                adjust_covariates=adjust_covariates,
            )
            rho, pval = stats.spearmanr(x, y)
            lr = stats.linregress(x, y)
            if not np.isfinite(rho) or not np.isfinite(pval):
                continue
            ci_stats = _bootstrap_rho_slope_ci(
                x=x,
                y=y,
                n_boot=bootstrap_iterations,
                ci=bootstrap_ci,
                seed=seed + (pair_idx * 97),
            )
            pair_idx += 1
            rows.append(
                {
                    "cell_type": cell_type,
                    "signature": sc_col.replace("score__", ""),
                    "n_donors": n_donors,
                    "age_min": float(t["age"].min()),
                    "age_max": float(t["age"].max()),
                    "age_span": age_span,
                    "spearman_rho": float(rho),
                    "pvalue": float(pval),
                    "spearman_rho_ci_low": ci_stats["spearman_rho_ci_low"],
                    "spearman_rho_ci_high": ci_stats["spearman_rho_ci_high"],
                    "slope_per_year": float(lr.slope),
                    "slope_per_year_ci_low": ci_stats["slope_per_year_ci_low"],
                    "slope_per_year_ci_high": ci_stats["slope_per_year_ci_high"],
                    "effect_per_10y": float(lr.slope * 10.0),
                    "effect_per_10y_ci_low": float(ci_stats["slope_per_year_ci_low"] * 10.0)
                    if np.isfinite(ci_stats["slope_per_year_ci_low"])
                    else np.nan,
                    "effect_per_10y_ci_high": float(ci_stats["slope_per_year_ci_high"] * 10.0)
                    if np.isfinite(ci_stats["slope_per_year_ci_high"])
                    else np.nan,
                    "slope_pvalue": float(lr.pvalue),
                    "r_squared": float(lr.rvalue**2),
                    "adjusted": bool(used_adjustment),
                    "adjusted_for": adjusted_for if used_adjustment else "",
                }
            )

    if not rows:
        return _empty_association_table()

    out = pd.DataFrame(rows)
    out["fdr"] = _bh_fdr(out["pvalue"].to_numpy())
    out["direction"] = np.where(out["effect_per_10y"] > 0, "increase_with_age", "decrease_with_age")
    out["fdr_significant"] = out["fdr"] < 0.05
    out = out.sort_values(
        ["fdr", "pvalue", "spearman_rho"], ascending=[True, True, False]
    ).reset_index(drop=True)
    out["rank"] = np.arange(1, out.shape[0] + 1)
    return out


def _plot_heatmap(assoc: pd.DataFrame, path: Path, dpi: int) -> None:
    if assoc.empty:
        _save_placeholder(
            path, "Signature-Age Associations", "No valid donor-level tests were available.", dpi
        )
        return

    pivot_rho = assoc.pivot_table(
        index="cell_type",
        columns="signature",
        values="spearman_rho",
        aggfunc="mean",
    ).fillna(0.0)
    pivot_fdr = assoc.pivot_table(
        index="cell_type",
        columns="signature",
        values="fdr",
        aggfunc="min",
    )
    if pivot_rho.empty:
        _save_placeholder(path, "Signature-Age Associations", "Association matrix is empty.", dpi)
        return

    sig_count = (pivot_fdr < 0.05).sum(axis=1).fillna(0)
    row_order = sig_count.sort_values(ascending=False).index
    pivot_rho = pivot_rho.loc[row_order]
    pivot_fdr = pivot_fdr.loc[row_order]

    fig_h = max(4, 0.35 * len(pivot_rho.index))
    fig_w = max(6, 1.2 * len(pivot_rho.columns))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(pivot_rho.to_numpy(), aspect="auto", cmap=DIVERGING_CMAP, vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(pivot_rho.columns)))
    ax.set_xticklabels(
        [_format_signature_label(c) for c in pivot_rho.columns], rotation=30, ha="right"
    )
    ax.set_yticks(np.arange(len(pivot_rho.index)))
    ax.set_yticklabels(pivot_rho.index)
    ax.set_title("Signature-age associations by immune cell population", pad=12)

    for i in range(pivot_rho.shape[0]):
        for j in range(pivot_rho.shape[1]):
            fdr = pivot_fdr.iloc[i, j]
            if pd.notna(fdr) and fdr < 0.05:
                ax.text(
                    j,
                    i,
                    "*",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=10,
                    fontweight="bold",
                )

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Spearman rho")
    fig.tight_layout()
    fdr_note = cbar.ax.text(
        1.15,
        -0.28,
        "* FDR < 0.05",
        transform=cbar.ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#4A4A4A",
        clip_on=False,
    )
    fig.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
        bbox_extra_artists=(fdr_note,),
        pad_inches=0.2,
        facecolor="white",
    )
    plt.close(fig)


def _plot_top_panels(
    agg: pd.DataFrame,
    assoc: pd.DataFrame,
    path: Path,
    top_n: int,
    min_donors_for_panel: int,
    significant_only: bool,
    dpi: int,
) -> None:
    if assoc.empty or agg.empty:
        _save_placeholder(
            path, "Top Signature-Age Associations", "No associations to display.", dpi
        )
        return

    candidates = _filter_assoc_for_panels(
        assoc=assoc,
        significant_only=significant_only,
        min_donors_for_panel=min_donors_for_panel,
    )
    if candidates.empty:
        filters = []
        if significant_only:
            filters.append("FDR<0.05")
        if min_donors_for_panel > 0:
            filters.append(f"n_donors>={min_donors_for_panel}")
        suffix = f" ({', '.join(filters)})" if filters else ""
        _save_placeholder(
            path,
            "Top Signature-Age Associations",
            f"No associations satisfy panel filters{suffix}.",
            dpi,
        )
        return

    ranked = _sort_assoc_for_panels(candidates).head(top_n)
    n = ranked.shape[0]
    n_cols = 2
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows), squeeze=False)
    axes_flat = axes.reshape(-1)

    for i, row in enumerate(ranked.itertuples(index=False)):
        ax = axes_flat[i]
        sig_col = f"score__{row.signature}"
        subset = agg[(agg["cell_type"] == row.cell_type)][["age", sig_col]].dropna()
        x = subset["age"].to_numpy(dtype=float)
        y = subset[sig_col].to_numpy(dtype=float)
        ax.scatter(x, y, s=20, alpha=0.75, color=PALETTE["primary"])
        if subset.shape[0] >= 2:
            lr = stats.linregress(x, y)
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            ys = lr.intercept + lr.slope * xs
            line_color = PALETTE["secondary"] if lr.slope >= 0 else PALETTE["danger"]
            ax.plot(xs, ys, linewidth=2, color=line_color)
        sig_label = _format_signature_label(str(row.signature))
        adjusted_note = ""
        if bool(getattr(row, "adjusted", False)):
            adjusted_for = str(getattr(row, "adjusted_for", "")).strip()
            adjusted_note = f", adj={adjusted_for}" if adjusted_for else ", adj"
        ax.set_title(
            f"{row.cell_type} | {sig_label}\n"
            f"rho={row.spearman_rho:.2f}, FDR={row.fdr:.2e}, effect/10y={row.effect_per_10y:.3f}{adjusted_note}"
        )
        ax.set_xlabel("Age")
        ax.set_ylabel("Mean signature score")
        style_axis(ax, grid="y")

    for ax in axes_flat[n:]:
        ax.axis("off")

    finalize_and_save(fig, path, dpi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--fig-heatmap", required=True)
    ap.add_argument("--fig-top", required=True)
    ap.add_argument("--table-scores", required=True)
    ap.add_argument("--table-assoc", required=True)
    ap.add_argument("--table-signature-meta", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    scfg = cfg.get("signature_age", {})
    signatures = _extract_signatures(cfg)

    dpi = int(scfg.get("dpi", cfg.get("report", {}).get("dpi", 150)))
    min_genes_present = int(scfg.get("min_genes_present", 3))
    min_cells_per_group = int(scfg.get("min_cells_per_group", 50))
    min_donors_per_celltype = int(scfg.get("min_donors_per_celltype", 8))
    min_age_span = float(scfg.get("min_age_span", 10.0))
    min_donors_for_panel = int(scfg.get("min_donors_for_panel", max(20, min_donors_per_celltype)))
    max_cells_for_scoring = int(scfg.get("max_cells_for_scoring", 200000))
    top_n_panels = int(scfg.get("top_n_panels", 8))
    plot_only_fdr_significant = bool(scfg.get("plot_only_fdr_significant", True))
    adjust_covariates = bool(scfg.get("adjust_covariates", True))
    bootstrap_iterations = int(
        scfg.get(
            "bootstrap_iterations", cfg.get("age_prediction", {}).get("bootstrap_iterations", 2000)
        )
    )
    bootstrap_ci = float(
        scfg.get("bootstrap_ci", cfg.get("age_prediction", {}).get("bootstrap_ci", 0.95))
    )
    sampling_strategy = str(scfg.get("sampling_strategy", "stratified")).lower()
    seed = int(cfg.get("run", {}).get("seed", 42))
    assume_log1p = bool(scfg.get("assume_log1p", False))
    target_sum = float(scfg.get("normalize_target_sum", 1e4))
    materialize_chunk_size = int(scfg.get("materialize_chunk_size", 5000))
    bootstrap_iterations = max(0, bootstrap_iterations)
    bootstrap_ci = float(np.clip(bootstrap_ci, 0.5, 0.999))

    ensure_dir(Path(args.fig_heatmap).parent)
    ensure_dir(Path(args.fig_top).parent)
    ensure_dir(Path(args.table_scores).parent)
    ensure_dir(Path(args.table_assoc).parent)
    if args.table_signature_meta:
        ensure_dir(Path(args.table_signature_meta).parent)
    apply_publication_style(dpi=dpi)

    adata_backed = ad.read_h5ad(args.inp, backed="r")
    attach_biological_replicates(adata_backed, cfg)
    colmap = _pick_columns(list(adata_backed.obs.columns), cfg)
    covariate_cols = _resolve_covariate_cols(list(adata_backed.obs.columns), cfg, colmap)
    missing = [k for k, v in colmap.items() if k in {"age", "donor", "cell_type"} and v is None]
    if missing:
        if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
            adata_backed.file.close()
        msg = f"Missing required columns: {', '.join(missing)}"
        require_placeholder_permission(cfg, msg)
        pd.DataFrame().to_csv(args.table_scores, index=False)
        _empty_association_table().to_csv(args.table_assoc, index=False)
        if args.table_signature_meta:
            pd.DataFrame(
                columns=[
                    "signature",
                    "n_genes_requested",
                    "n_genes_present",
                    "present_genes",
                    "missing_genes",
                ]
            ).to_csv(
                args.table_signature_meta,
                index=False,
            )
        _save_placeholder(Path(args.fig_heatmap), "Signature-Age Associations", msg, dpi)
        _save_placeholder(Path(args.fig_top), "Top Signature-Age Associations", msg, dpi)
        return

    sel_idx = _subset_cells(
        adata_backed,
        age_col=str(colmap["age"]),
        donor_col=str(colmap["donor"]),
        celltype_col=str(colmap["cell_type"]),
        max_cells=max_cells_for_scoring,
        min_cells_per_group=min_cells_per_group,
        seed=seed,
        sampling_strategy=sampling_strategy,
    )
    print(
        f"[signature_age] selected {len(sel_idx):,} cells for scoring",
        flush=True,
    )
    if len(sel_idx) == 0:
        if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
            adata_backed.file.close()
        require_placeholder_permission(
            cfg, "Signature analysis has no eligible cells after filtering."
        )
        pd.DataFrame().to_csv(args.table_scores, index=False)
        _empty_association_table().to_csv(args.table_assoc, index=False)
        if args.table_signature_meta:
            pd.DataFrame(
                columns=[
                    "signature",
                    "n_genes_requested",
                    "n_genes_present",
                    "present_genes",
                    "missing_genes",
                ]
            ).to_csv(
                args.table_signature_meta,
                index=False,
            )
        _save_placeholder(
            Path(args.fig_heatmap),
            "Signature-Age Associations",
            "No eligible cells after filtering.",
            dpi,
        )
        _save_placeholder(
            Path(args.fig_top),
            "Top Signature-Age Associations",
            "No eligible cells after filtering.",
            dpi,
        )
        return

    adata = _materialize_signature_genes(
        adata_backed,
        selected_indices=sel_idx,
        signatures=signatures,
        assume_log1p=assume_log1p,
        target_sum=target_sum,
        chunk_size=materialize_chunk_size,
    )
    if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
        adata_backed.file.close()

    print("[signature_age] computing signature scores", flush=True)
    score_df, signature_meta = _mean_expression_per_signature(
        adata,
        signatures=signatures,
        min_genes_present=min_genes_present,
    )
    score_cols = [c for c in score_df.columns if c.startswith("score__")]

    obs = adata.obs.copy()
    obs["age"] = pd.to_numeric(obs[str(colmap["age"])], errors="coerce")
    obs["donor_id"] = obs[str(colmap["donor"])].astype("string").str.strip()
    obs["cell_type"] = obs[str(colmap["cell_type"])].astype("string").str.strip()
    obs["donor_id"] = obs["donor_id"].mask(obs["donor_id"].eq(""))
    obs["cell_type"] = obs["cell_type"].mask(obs["cell_type"].eq(""))
    for cov in covariate_cols:
        if cov in adata.obs.columns:
            obs[cov] = adata.obs[cov]

    base_cols = ["donor_id", "age", "cell_type"] + [c for c in covariate_cols if c in obs.columns]
    obs = pd.concat([obs[base_cols], score_df], axis=1)
    obs = obs.dropna(subset=["age", "donor_id", "cell_type"])
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    validate_replicate_covariates(
        obs,
        replicate_col="donor_id",
        covariate_cols=covariate_cols,
        context="signature_age",
    )

    agg = _aggregate_scores(
        obs,
        score_cols=score_cols,
        min_cells_per_group=min_cells_per_group,
        covariate_cols=covariate_cols,
    )
    assoc = _association_table(
        agg,
        score_cols=score_cols,
        min_donors_per_celltype=min_donors_per_celltype,
        min_age_span=min_age_span,
        covariate_cols=covariate_cols,
        adjust_covariates=adjust_covariates,
        bootstrap_iterations=bootstrap_iterations,
        bootstrap_ci=bootstrap_ci,
        seed=seed,
    )
    print(
        "[signature_age] "
        f"computed {len(assoc):,} association tests from "
        f"{agg['donor_id'].nunique():,} replicates",
        flush=True,
    )

    agg = agg.sort_values(["cell_type", "donor_id"]).reset_index(drop=True)
    grouping_col = str(scfg.get("donor_col", "donor_id"))
    agg_out = agg.rename(columns={"donor_id": grouping_col})
    assoc["grouping_id_column"] = grouping_col
    agg_out.to_csv(args.table_scores, index=False)
    assoc.to_csv(args.table_assoc, index=False)
    if args.table_signature_meta:
        signature_meta.to_csv(args.table_signature_meta, index=False)

    _plot_heatmap(assoc, path=Path(args.fig_heatmap), dpi=dpi)
    _plot_top_panels(
        agg,
        assoc,
        path=Path(args.fig_top),
        top_n=top_n_panels,
        min_donors_for_panel=min_donors_for_panel,
        significant_only=plot_only_fdr_significant,
        dpi=dpi,
    )


if __name__ == "__main__":
    main()
