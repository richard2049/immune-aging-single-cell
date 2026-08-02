from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .plot_style import DIVERGING_CMAP, PALETTE, apply_publication_style, finalize_and_save, save_placeholder, style_axis
from .utils import ensure_dir, load_config


def _safe_bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"})


def _safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def _format_signature_label(name: str) -> str:
    mapping = {
        "ifn_response": "IFN response",
        "inflammatory_nfkb": "Inflammatory (NF-kB)",
        "mitochondrial_stress": "Mitochondrial stress",
        "proteostasis_upr": "Proteostasis (UPR)",
        "sasp_proxy": "SASP proxy",
    }
    if name in mapping:
        return mapping[name]
    return str(name).replace("_", " ").strip().title()


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _prepare_comp_effects(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "cell_type" not in df.columns:
        return pd.DataFrame(columns=["signal_id", "label", "effect", "ci_low", "ci_high", "fdr", "significant"])
    out = pd.DataFrame()
    out["signal_id"] = df["cell_type"].astype(str)
    out["label"] = df["cell_type"].astype(str)
    out["effect"] = pd.to_numeric(df.get("slope_per_10y"), errors="coerce")
    out["ci_low"] = pd.to_numeric(df.get("slope_per_10y_ci_low"), errors="coerce")
    out["ci_high"] = pd.to_numeric(df.get("slope_per_10y_ci_high"), errors="coerce")
    out["fdr"] = pd.to_numeric(df.get("spearman_fdr"), errors="coerce")
    if "fdr_significant" in df.columns:
        out["significant"] = _safe_bool_series(df["fdr_significant"])
    else:
        out["significant"] = out["fdr"] < 0.05
    return out.dropna(subset=["effect"]).drop_duplicates(subset=["signal_id"])


def _prepare_sig_effects(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or any(c not in df.columns for c in ["cell_type", "signature"]):
        return pd.DataFrame(columns=["signal_id", "label", "effect", "ci_low", "ci_high", "fdr", "significant"])
    out = pd.DataFrame()
    out["cell_type"] = df["cell_type"].astype(str)
    out["signature"] = df["signature"].astype(str)
    out["signal_id"] = out["cell_type"] + "||" + out["signature"]
    out["label"] = out["cell_type"] + " | " + out["signature"].map(_format_signature_label)
    out["effect"] = pd.to_numeric(df.get("effect_per_10y"), errors="coerce")
    out["ci_low"] = pd.to_numeric(df.get("effect_per_10y_ci_low"), errors="coerce")
    out["ci_high"] = pd.to_numeric(df.get("effect_per_10y_ci_high"), errors="coerce")
    out["fdr"] = pd.to_numeric(df.get("fdr"), errors="coerce")
    if "fdr_significant" in df.columns:
        out["significant"] = _safe_bool_series(df["fdr_significant"])
    else:
        out["significant"] = out["fdr"] < 0.05
    return out.dropna(subset=["effect"]).drop_duplicates(subset=["signal_id"])


def _select_top_signals(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    if df.empty:
        return df
    ranked = df.copy()
    ranked["fdr_rank"] = ranked["fdr"].fillna(np.inf)
    ranked["abs_effect"] = ranked["effect"].abs()
    ranked = ranked.sort_values(["significant", "fdr_rank", "abs_effect"], ascending=[False, True, False]).reset_index(drop=True)
    return ranked.head(int(max(top_n, 1))).reset_index(drop=True)


def _plot_forest_panel(ax, df: pd.DataFrame, panel_title: str, xlabel: str) -> None:
    if df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No signals available", ha="center", va="center")
        ax.set_title(panel_title)
        return

    plot_df = df.copy().iloc[::-1].reset_index(drop=True)
    y = np.arange(plot_df.shape[0])
    ax.axvline(0.0, color=PALETTE["muted"], linewidth=1.2, linestyle="--")
    for i, row in plot_df.iterrows():
        effect = float(row["effect"])
        lo = float(row["ci_low"]) if np.isfinite(row["ci_low"]) else np.nan
        hi = float(row["ci_high"]) if np.isfinite(row["ci_high"]) else np.nan
        color = PALETTE["secondary"] if effect >= 0 else PALETTE["danger"]
        if np.isfinite(lo) and np.isfinite(hi):
            xerr = np.array([[max(effect - lo, 0.0)], [max(hi - effect, 0.0)]])
            ax.errorbar(effect, i, xerr=xerr, fmt="o", color=color, ecolor=color, capsize=3, markersize=5, linewidth=1.5)
        else:
            ax.scatter(effect, i, s=24, color=color)
        if bool(row.get("significant", False)):
            ax.scatter(effect, i, s=70, facecolors="none", edgecolors="#222222", linewidths=1.0)

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["label"].tolist())
    ax.set_title(panel_title)
    ax.set_xlabel(xlabel)
    style_axis(ax, grid="x")


def _plot_effect_ci_forest(
    comp_top: pd.DataFrame,
    sig_top: pd.DataFrame,
    path: Path,
    dpi: int,
) -> None:
    if comp_top.empty and sig_top.empty:
        save_placeholder(path, "Supplementary: Effect Sizes with 95% CIs", "No composition/signature effects available.", dpi)
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.3]})
    _plot_forest_panel(
        axes[0],
        comp_top,
        panel_title="Cell-type composition trends",
        xlabel="Effect per 10y (fraction units)",
    )
    _plot_forest_panel(
        axes[1],
        sig_top,
        panel_title="Signature-age associations",
        xlabel="Effect per 10y (signature-score units)",
    )
    fig.suptitle("Supplementary: Donor-level age effects with 95% bootstrap CIs", y=1.02)
    finalize_and_save(fig, path, dpi)


def _scenario_label(row: pd.Series) -> str:
    if {"adjust_covariates", "drop_sparse_age_bins", "min_cells_per_group"}.issubset(set(row.index)):
        a = int(_as_bool(row["adjust_covariates"]))
        d = int(_as_bool(row["drop_sparse_age_bins"]))
        try:
            m = int(row["min_cells_per_group"])
        except Exception:
            m = -1
        return f"A{a}|D{d}|M{m}"
    return str(row.get("scenario", "scenario"))


def _read_effect_lookup(path: Path, mode: str) -> pd.DataFrame:
    df = _safe_read_csv(path)
    if df.empty:
        return pd.DataFrame(columns=["signal_id", "effect", "significant"])
    if mode == "composition":
        out = _prepare_comp_effects(df)
    else:
        out = _prepare_sig_effects(df)
    keep = out[["signal_id", "effect", "significant"]].copy()
    return keep.drop_duplicates(subset=["signal_id"])


def _build_stability_matrix(
    manifest: pd.DataFrame,
    selected: pd.DataFrame,
    mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if manifest.empty or selected.empty:
        return pd.DataFrame(), pd.DataFrame()

    label_by_signal = selected.set_index("signal_id")["label"].to_dict()
    scenario_labels = manifest.apply(_scenario_label, axis=1).tolist()
    values = pd.DataFrame(np.nan, index=selected["signal_id"].tolist(), columns=scenario_labels, dtype=float)
    signif = pd.DataFrame(False, index=selected["signal_id"].tolist(), columns=scenario_labels, dtype=bool)

    for _, row in manifest.iterrows():
        col = _scenario_label(row)
        table_col = "composition_trends_table" if mode == "composition" else "signature_assoc_table"
        table_path = Path(str(row.get(table_col, "")))
        lookup = _read_effect_lookup(table_path, mode=mode)
        if lookup.empty:
            continue
        lookup = lookup.set_index("signal_id")
        for signal_id in values.index:
            if signal_id not in lookup.index:
                continue
            values.loc[signal_id, col] = float(lookup.loc[signal_id, "effect"])
            signif.loc[signal_id, col] = bool(lookup.loc[signal_id, "significant"])

    values.index = [label_by_signal.get(k, k) for k in values.index]
    signif.index = values.index
    return values, signif


def _plot_stability_panel(ax, values: pd.DataFrame, signif: pd.DataFrame, title: str) -> None:
    if values.empty or values.shape[1] < 2:
        ax.axis("off")
        ax.text(0.5, 0.5, "Run sensitivity scenarios to populate this panel", ha="center", va="center")
        ax.set_title(title)
        return

    arr = values.to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
    vmax = max(vmax, 1e-6)
    masked = np.ma.masked_invalid(arr)
    cmap = plt.get_cmap(DIVERGING_CMAP).copy()
    cmap.set_bad(color="#EAEAEA")
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(values.shape[1]))
    ax.set_xticklabels(values.columns.tolist(), rotation=35, ha="right")
    ax.set_yticks(np.arange(values.shape[0]))
    ax.set_yticklabels(values.index.tolist())
    ax.set_title(title)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if bool(signif.iloc[i, j]):
                ax.text(j, i, "*", ha="center", va="center", color="black", fontsize=9, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.75)
    cbar.set_label("Effect per 10y")


def _plot_stability_heatmap(
    manifest_path: Path,
    comp_top: pd.DataFrame,
    sig_top: pd.DataFrame,
    path: Path,
    dpi: int,
) -> None:
    manifest = _safe_read_csv(manifest_path)
    if manifest.empty:
        save_placeholder(
            path,
            "Supplementary: Sensitivity Stability",
            "Missing sensitivity manifest. Run sensitivity_age first.",
            dpi,
        )
        return

    if {"source", "adjust_covariates", "drop_sparse_age_bins", "min_cells_per_group"}.issubset(set(manifest.columns)):
        source_order = {"baseline": 0, "adjust_covariates": 1, "drop_sparse_age_bins": 2, "min_cells_per_group": 3}
        manifest = manifest.copy()
        manifest["_src_order"] = manifest["source"].map(lambda x: source_order.get(str(x), 9))
        manifest = manifest.sort_values(
            ["_src_order", "adjust_covariates", "drop_sparse_age_bins", "min_cells_per_group", "scenario"]
        ).drop(columns="_src_order")
    elif "scenario" in manifest.columns:
        manifest = manifest.sort_values("scenario")

    comp_vals, comp_sig = _build_stability_matrix(manifest=manifest, selected=comp_top, mode="composition")
    sig_vals, sig_sig = _build_stability_matrix(manifest=manifest, selected=sig_top, mode="signature")

    fig, axes = plt.subplots(2, 1, figsize=(13, 8), gridspec_kw={"height_ratios": [1, 1.2]})
    _plot_stability_panel(
        axes[0],
        comp_vals,
        comp_sig,
        title="Composition effect stability across sensitivity scenarios",
    )
    _plot_stability_panel(
        axes[1],
        sig_vals,
        sig_sig,
        title="Signature effect stability across sensitivity scenarios",
    )
    fig.text(
        0.01,
        0.01,
        "Scenario code: A=adjust_covariates, D=drop_sparse_age_bins, M=min_cells_per_group. Asterisks: FDR<0.05.",
        fontsize=9,
        color="#333333",
    )
    finalize_and_save(fig, path, dpi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--table-composition", required=True)
    ap.add_argument("--table-signature", required=True)
    ap.add_argument("--sensitivity-manifest", required=False, default="")
    ap.add_argument("--fig-effect-ci", required=True)
    ap.add_argument("--fig-stability", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    scfg = cfg.get("supplementary_age", {})
    dpi = int(scfg.get("dpi", cfg.get("report", {}).get("dpi", 180)))
    top_n_comp = int(scfg.get("top_n_composition", 8))
    top_n_sig = int(scfg.get("top_n_signature", 8))

    ensure_dir(Path(args.fig_effect_ci).parent)
    ensure_dir(Path(args.fig_stability).parent)
    apply_publication_style(dpi=dpi)

    comp_df = _safe_read_csv(Path(args.table_composition))
    sig_df = _safe_read_csv(Path(args.table_signature))

    comp_effects = _prepare_comp_effects(comp_df)
    sig_effects = _prepare_sig_effects(sig_df)
    comp_top = _select_top_signals(comp_effects, top_n=top_n_comp)
    sig_top = _select_top_signals(sig_effects, top_n=top_n_sig)

    _plot_effect_ci_forest(
        comp_top=comp_top,
        sig_top=sig_top,
        path=Path(args.fig_effect_ci),
        dpi=dpi,
    )
    _plot_stability_heatmap(
        manifest_path=Path(args.sensitivity_manifest) if args.sensitivity_manifest else Path("__missing__"),
        comp_top=comp_top,
        sig_top=sig_top,
        path=Path(args.fig_stability),
        dpi=dpi,
    )


if __name__ == "__main__":
    main()
