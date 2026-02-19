from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

from .plot_style import DIVERGING_CMAP, PALETTE, apply_publication_style, finalize_and_save, save_placeholder, style_axis
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
        donor = _first_existing(obs_cols, ["donor_id", "donor_id_y", "donor_id_x", "donor", "sample_id"])
    if celltype not in obs_cols:
        celltype = _first_existing(obs_cols, ["cell_type", "majority_voting", "predicted_labels", "leiden"])
    if sex not in obs_cols:
        sex = _first_existing(obs_cols, ["sex"])

    return {"age": age, "donor": donor, "cell_type": celltype, "sex": sex}


def _subset_cells(adata, age_col: str, donor_col: str, celltype_col: str, max_cells: int, seed: int) -> np.ndarray:
    obs = adata.obs
    keep = (
        pd.to_numeric(obs[age_col], errors="coerce").notna()
        & obs[donor_col].notna()
        & obs[celltype_col].notna()
    )
    idx = np.flatnonzero(keep.to_numpy())
    if len(idx) <= max_cells:
        return idx

    rng = np.random.default_rng(seed)
    chosen = rng.choice(idx, size=max_cells, replace=False)
    chosen.sort()
    return chosen


def _mean_expression_per_signature(adata, signatures: dict[str, list[str]], min_genes_present: int) -> tuple[pd.DataFrame, pd.DataFrame]:
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


def _aggregate_scores(obs: pd.DataFrame, score_cols: list[str], min_cells_per_group: int) -> pd.DataFrame:
    grouped = (
        obs.groupby(["donor_id", "age", "sex", "cell_type"], observed=False)
        .agg(n_cells=("cell_type", "size"), **{c: (c, "mean") for c in score_cols})
        .reset_index()
    )
    grouped = grouped[grouped["n_cells"] >= min_cells_per_group].copy()
    return grouped


def _association_table(
    agg: pd.DataFrame,
    score_cols: list[str],
    min_donors_per_celltype: int,
    min_age_span: float,
) -> pd.DataFrame:
    rows = []
    for cell_type, gct in agg.groupby("cell_type", observed=False):
        for sc_col in score_cols:
            t = gct[["age", sc_col]].dropna().copy()
            if t.shape[0] < min_donors_per_celltype:
                continue
            age_span = float(t["age"].max() - t["age"].min())
            if age_span < min_age_span:
                continue
            rho, pval = stats.spearmanr(t["age"].to_numpy(), t[sc_col].to_numpy())
            lr = stats.linregress(t["age"].to_numpy(), t[sc_col].to_numpy())
            rows.append(
                {
                    "cell_type": cell_type,
                    "signature": sc_col.replace("score__", ""),
                    "n_donors": int(t.shape[0]),
                    "age_min": float(t["age"].min()),
                    "age_max": float(t["age"].max()),
                    "age_span": age_span,
                    "spearman_rho": float(rho),
                    "pvalue": float(pval),
                    "slope_per_year": float(lr.slope),
                    "effect_per_10y": float(lr.slope * 10.0),
                    "slope_pvalue": float(lr.pvalue),
                    "r_squared": float(lr.rvalue**2),
                }
            )

    if not rows:
        cols = [
            "cell_type",
            "signature",
            "n_donors",
            "age_min",
            "age_max",
            "age_span",
            "spearman_rho",
            "pvalue",
            "slope_per_year",
            "effect_per_10y",
            "slope_pvalue",
            "r_squared",
            "fdr",
            "direction",
            "fdr_significant",
            "rank",
        ]
        return pd.DataFrame(columns=cols)

    out = pd.DataFrame(rows)
    out["fdr"] = _bh_fdr(out["pvalue"].to_numpy())
    out["direction"] = np.where(out["effect_per_10y"] > 0, "increase_with_age", "decrease_with_age")
    out["fdr_significant"] = out["fdr"] < 0.05
    out = out.sort_values(["fdr", "pvalue", "spearman_rho"], ascending=[True, True, False]).reset_index(drop=True)
    out["rank"] = np.arange(1, out.shape[0] + 1)
    return out


def _plot_heatmap(assoc: pd.DataFrame, path: Path, dpi: int) -> None:
    if assoc.empty:
        _save_placeholder(path, "Signature-Age Associations", "No valid donor-level tests were available.", dpi)
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
    ax.set_xticklabels(pivot_rho.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pivot_rho.index)))
    ax.set_yticklabels(pivot_rho.index)
    ax.set_title("Signature-age associations (Spearman rho)\n'*' marks FDR < 0.05")

    for i in range(pivot_rho.shape[0]):
        for j in range(pivot_rho.shape[1]):
            fdr = pivot_fdr.iloc[i, j]
            if pd.notna(fdr) and fdr < 0.05:
                ax.text(j, i, "*", ha="center", va="center", color="black", fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Spearman rho")
    finalize_and_save(fig, path, dpi)


def _plot_top_panels(agg: pd.DataFrame, assoc: pd.DataFrame, path: Path, top_n: int, dpi: int) -> None:
    if assoc.empty or agg.empty:
        _save_placeholder(path, "Top Signature-Age Associations", "No associations to display.", dpi)
        return

    ranked = assoc.sort_values(["fdr", "pvalue", "spearman_rho"], ascending=[True, True, False]).head(top_n)
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
        ax.set_title(
            f"{row.cell_type} | {row.signature}\n"
            f"rho={row.spearman_rho:.2f}, FDR={row.fdr:.2e}, effect/10y={row.effect_per_10y:.3f}"
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
    max_cells_for_scoring = int(scfg.get("max_cells_for_scoring", 200000))
    top_n_panels = int(scfg.get("top_n_panels", 8))
    seed = int(cfg.get("run", {}).get("seed", 42))
    assume_log1p = bool(scfg.get("assume_log1p", False))
    target_sum = float(scfg.get("normalize_target_sum", 1e4))

    ensure_dir(Path(args.fig_heatmap).parent)
    ensure_dir(Path(args.fig_top).parent)
    ensure_dir(Path(args.table_scores).parent)
    ensure_dir(Path(args.table_assoc).parent)
    if args.table_signature_meta:
        ensure_dir(Path(args.table_signature_meta).parent)
    apply_publication_style(dpi=dpi)

    adata_backed = sc.read_h5ad(args.inp, backed="r")
    colmap = _pick_columns(list(adata_backed.obs.columns), cfg)
    missing = [k for k, v in colmap.items() if k in {"age", "donor", "cell_type"} and v is None]
    if missing:
        if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
            adata_backed.file.close()
        msg = f"Missing required columns: {', '.join(missing)}"
        pd.DataFrame().to_csv(args.table_scores, index=False)
        pd.DataFrame().to_csv(args.table_assoc, index=False)
        if args.table_signature_meta:
            pd.DataFrame(columns=["signature", "n_genes_requested", "n_genes_present", "present_genes", "missing_genes"]).to_csv(
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
        seed=seed,
    )
    if len(sel_idx) == 0:
        if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
            adata_backed.file.close()
        pd.DataFrame().to_csv(args.table_scores, index=False)
        pd.DataFrame().to_csv(args.table_assoc, index=False)
        if args.table_signature_meta:
            pd.DataFrame(columns=["signature", "n_genes_requested", "n_genes_present", "present_genes", "missing_genes"]).to_csv(
                args.table_signature_meta,
                index=False,
            )
        _save_placeholder(Path(args.fig_heatmap), "Signature-Age Associations", "No eligible cells after filtering.", dpi)
        _save_placeholder(Path(args.fig_top), "Top Signature-Age Associations", "No eligible cells after filtering.", dpi)
        return

    adata = adata_backed[sel_idx].to_memory()
    if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
        adata_backed.file.close()
    if not assume_log1p:
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)

    score_df, signature_meta = _mean_expression_per_signature(
        adata,
        signatures=signatures,
        min_genes_present=min_genes_present,
    )
    score_cols = [c for c in score_df.columns if c.startswith("score__")]

    obs = adata.obs.copy()
    obs["age"] = pd.to_numeric(obs[str(colmap["age"])], errors="coerce")
    obs["donor_id"] = obs[str(colmap["donor"])].astype(str)
    obs["cell_type"] = obs[str(colmap["cell_type"])].astype(str)
    if colmap["sex"] is not None:
        obs["sex"] = obs[str(colmap["sex"])].astype(str)
    else:
        obs["sex"] = "NA"

    obs = pd.concat([obs[["donor_id", "age", "sex", "cell_type"]], score_df], axis=1)
    obs = obs.dropna(subset=["age", "donor_id", "cell_type"])

    agg = _aggregate_scores(obs, score_cols=score_cols, min_cells_per_group=min_cells_per_group)
    assoc = _association_table(
        agg,
        score_cols=score_cols,
        min_donors_per_celltype=min_donors_per_celltype,
        min_age_span=min_age_span,
    )

    agg = agg.sort_values(["cell_type", "donor_id"]).reset_index(drop=True)
    agg.to_csv(args.table_scores, index=False)
    assoc.to_csv(args.table_assoc, index=False)
    if args.table_signature_meta:
        signature_meta.to_csv(args.table_signature_meta, index=False)

    _plot_heatmap(assoc, path=Path(args.fig_heatmap), dpi=dpi)
    _plot_top_panels(agg, assoc, path=Path(args.fig_top), top_n=top_n_panels, dpi=dpi)


if __name__ == "__main__":
    main()
