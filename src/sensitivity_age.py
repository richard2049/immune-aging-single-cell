from __future__ import annotations

import argparse
import copy
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .utils import ensure_dir, load_config


def _as_bool_list(values: object, default: list[bool]) -> list[bool]:
    if not isinstance(values, list) or not values:
        return default
    out: list[bool] = []
    for v in values:
        if isinstance(v, bool):
            val = v
        else:
            sval = str(v).strip().lower()
            val = sval in {"1", "true", "yes", "y", "on"}
        if val not in out:
            out.append(val)
    return out if out else default


def _as_int_list(values: object, default: list[int]) -> list[int]:
    if not isinstance(values, list) or not values:
        return default
    out: list[int] = []
    for v in values:
        try:
            val = int(v)
        except Exception:
            continue
        if val > 0 and val not in out:
            out.append(val)
    return out if out else default


def _scenario_name(
    adjust_covariates: bool,
    drop_sparse_age_bins: bool,
    min_cells_per_group: int,
) -> str:
    return (
        f"adj_{int(adjust_covariates)}"
        f"__drop_sparse_{int(drop_sparse_age_bins)}"
        f"__min_cells_{int(min_cells_per_group)}"
    )


def _build_scenarios(cfg: dict) -> list[dict]:
    ccfg = cfg.get("composition_age", {})
    scfg = cfg.get("signature_age", {})
    sens = cfg.get("sensitivity_age", {})

    base_adjust = bool(ccfg.get("adjust_covariates", scfg.get("adjust_covariates", True)))
    base_drop_sparse = bool(ccfg.get("drop_sparse_age_bins", True))
    base_min_cells = int(scfg.get("min_cells_per_group", 120))

    adjust_values = _as_bool_list(sens.get("adjust_covariates_values"), [True, False])
    drop_values = _as_bool_list(sens.get("drop_sparse_age_bins_values"), [True, False])
    min_cells_values = _as_int_list(sens.get("min_cells_per_group_values"), [80, 120, 160])
    include_baseline = bool(sens.get("include_baseline", True))

    scenarios: list[dict] = []
    if include_baseline:
        scenarios.append(
            {
                "scenario": _scenario_name(base_adjust, base_drop_sparse, base_min_cells),
                "adjust_covariates": base_adjust,
                "drop_sparse_age_bins": base_drop_sparse,
                "min_cells_per_group": base_min_cells,
                "source": "baseline",
            }
        )

    for val in adjust_values:
        scenarios.append(
            {
                "scenario": _scenario_name(bool(val), base_drop_sparse, base_min_cells),
                "adjust_covariates": bool(val),
                "drop_sparse_age_bins": base_drop_sparse,
                "min_cells_per_group": base_min_cells,
                "source": "adjust_covariates",
            }
        )

    for val in drop_values:
        scenarios.append(
            {
                "scenario": _scenario_name(base_adjust, bool(val), base_min_cells),
                "adjust_covariates": base_adjust,
                "drop_sparse_age_bins": bool(val),
                "min_cells_per_group": base_min_cells,
                "source": "drop_sparse_age_bins",
            }
        )

    for val in min_cells_values:
        scenarios.append(
            {
                "scenario": _scenario_name(base_adjust, base_drop_sparse, int(val)),
                "adjust_covariates": base_adjust,
                "drop_sparse_age_bins": base_drop_sparse,
                "min_cells_per_group": int(val),
                "source": "min_cells_per_group",
            }
        )

    deduped: list[dict] = []
    seen: set[tuple[bool, bool, int]] = set()
    for s in scenarios:
        key = (bool(s["adjust_covariates"]), bool(s["drop_sparse_age_bins"]), int(s["min_cells_per_group"]))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(s)
    return deduped


def _run_command(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _scenario_metrics(trends_path: Path, assoc_path: Path) -> dict[str, float]:
    metrics: dict[str, float] = {
        "n_trend_tests": np.nan,
        "n_trend_fdr_lt_0_05": np.nan,
        "median_abs_slope_per_10y": np.nan,
        "n_signature_tests": np.nan,
        "n_signature_fdr_lt_0_05": np.nan,
        "median_abs_signature_effect_per_10y": np.nan,
    }

    if trends_path.exists():
        trends = pd.read_csv(trends_path)
        metrics["n_trend_tests"] = float(trends.shape[0])
        if "fdr_significant" in trends.columns:
            sig = pd.Series(trends["fdr_significant"]).astype(str).str.lower().isin({"true", "1", "yes"})
            metrics["n_trend_fdr_lt_0_05"] = float(sig.sum())
        if "slope_per_10y" in trends.columns and not trends.empty:
            metrics["median_abs_slope_per_10y"] = float(pd.to_numeric(trends["slope_per_10y"], errors="coerce").abs().median())

    if assoc_path.exists():
        assoc = pd.read_csv(assoc_path)
        metrics["n_signature_tests"] = float(assoc.shape[0])
        if "fdr_significant" in assoc.columns:
            sig = pd.Series(assoc["fdr_significant"]).astype(str).str.lower().isin({"true", "1", "yes"})
            metrics["n_signature_fdr_lt_0_05"] = float(sig.sum())
        if "effect_per_10y" in assoc.columns and not assoc.empty:
            metrics["median_abs_signature_effect_per_10y"] = float(pd.to_numeric(assoc["effect_per_10y"], errors="coerce").abs().median())

    return metrics


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--done", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ensure_dir(Path(args.manifest).parent)
    ensure_dir(Path(args.summary).parent)
    ensure_dir(Path(args.done).parent)

    scenarios = _build_scenarios(cfg)
    sensitivity_bootstrap_iterations = int(
        cfg.get("sensitivity_age", {}).get("bootstrap_iterations", 500)
    )
    manifest_rows: list[dict] = []
    summary_rows: list[dict] = []

    for i, scenario in enumerate(scenarios):
        scenario_name = str(scenario["scenario"])
        print(
            "[sensitivity_age] "
            f"scenario {i + 1}/{len(scenarios)}: {scenario_name}",
            flush=True,
        )
        scenario_dir = outdir / scenario_name
        figdir = scenario_dir / "figures"
        tabledir = scenario_dir / "tables"
        ensure_dir(figdir)
        ensure_dir(tabledir)

        cfg_s = copy.deepcopy(cfg)
        cfg_s.setdefault("composition_age", {})
        cfg_s.setdefault("signature_age", {})
        cfg_s["composition_age"]["adjust_covariates"] = bool(scenario["adjust_covariates"])
        cfg_s["signature_age"]["adjust_covariates"] = bool(scenario["adjust_covariates"])
        cfg_s["composition_age"]["drop_sparse_age_bins"] = bool(scenario["drop_sparse_age_bins"])
        cfg_s["signature_age"]["min_cells_per_group"] = int(scenario["min_cells_per_group"])
        cfg_s["composition_age"]["bootstrap_iterations"] = (
            sensitivity_bootstrap_iterations
        )
        cfg_s["signature_age"]["bootstrap_iterations"] = (
            sensitivity_bootstrap_iterations
        )
        cfg_s.setdefault("sensitivity_age", {})
        cfg_s["sensitivity_age"]["active_scenario"] = scenario_name
        cfg_s["sensitivity_age"]["active_scenario_index"] = int(i)

        scenario_cfg = scenario_dir / "config.sensitivity.yml"
        with scenario_cfg.open("w", encoding="utf-8") as f:
            yaml.safe_dump(cfg_s, f, sort_keys=False)

        comp_fig = figdir / "age_celltype_composition_by_bin.png"
        comp_trends_fig = figdir / "age_celltype_top_trends.png"
        comp_tbl = tabledir / "age_celltype_fraction_by_donor.csv"
        comp_stats = tabledir / "age_celltype_trend_stats.csv"

        sig_heat = figdir / "signature_age_heatmap.png"
        sig_top = figdir / "signature_age_top_associations.png"
        sig_scores = tabledir / "signature_scores_by_donor_celltype.csv"
        sig_assoc = tabledir / "signature_age_associations.csv"
        sig_meta = tabledir / "signature_gene_coverage.csv"

        _run_command(
            [
                sys.executable,
                "-m",
                "src.composition_age",
                "--config",
                str(scenario_cfg),
                "--inp",
                str(args.inp),
                "--fig-fractions",
                str(comp_fig),
                "--fig-trends",
                str(comp_trends_fig),
                "--table-donor-fractions",
                str(comp_tbl),
                "--table-trends",
                str(comp_stats),
            ]
        )
        _run_command(
            [
                sys.executable,
                "-m",
                "src.signature_age",
                "--config",
                str(scenario_cfg),
                "--inp",
                str(args.inp),
                "--fig-heatmap",
                str(sig_heat),
                "--fig-top",
                str(sig_top),
                "--table-scores",
                str(sig_scores),
                "--table-assoc",
                str(sig_assoc),
                "--table-signature-meta",
                str(sig_meta),
            ]
        )

        manifest_rows.append(
            {
                "scenario": scenario_name,
                "source": str(scenario["source"]),
                "adjust_covariates": bool(scenario["adjust_covariates"]),
                "drop_sparse_age_bins": bool(scenario["drop_sparse_age_bins"]),
                "min_cells_per_group": int(scenario["min_cells_per_group"]),
                "scenario_dir": str(scenario_dir),
                "scenario_config": str(scenario_cfg),
                "composition_trends_table": str(comp_stats),
                "signature_assoc_table": str(sig_assoc),
            }
        )
        summary_rows.append(
            {
                "scenario": scenario_name,
                "adjust_covariates": bool(scenario["adjust_covariates"]),
                "drop_sparse_age_bins": bool(scenario["drop_sparse_age_bins"]),
                "min_cells_per_group": int(scenario["min_cells_per_group"]),
                **_scenario_metrics(comp_stats, sig_assoc),
            }
        )

    manifest = pd.DataFrame(manifest_rows).sort_values("scenario").reset_index(drop=True)
    summary = pd.DataFrame(summary_rows).sort_values("scenario").reset_index(drop=True)
    manifest.to_csv(args.manifest, index=False)
    summary.to_csv(args.summary, index=False)
    Path(args.done).write_text("sensitivity_age_complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()
