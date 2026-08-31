from __future__ import annotations

import argparse
import os
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

from .scvi_train import _index_sha256, _resolved_config_sha256
from .utils import load_config


def _required_clustering_arguments(
    clustering: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    groups = {
        "neighbors_args": {
            "n_neighbors",
            "n_pcs",
            "knn",
            "method",
            "transformer",
            "metric",
            "metric_kwds",
        },
        "umap_args": {
            "min_dist",
            "spread",
            "n_components",
            "maxiter",
            "alpha",
            "gamma",
            "negative_sample_rate",
            "init_pos",
            "a",
            "b",
            "method",
        },
        "leiden_args": {
            "directed",
            "use_weights",
            "n_iterations",
            "neighbors_key",
            "flavor",
        },
    }
    resolved: dict[str, dict[str, Any]] = {}
    missing: dict[str, list[str]] = {}
    for name, required in groups.items():
        values = dict(clustering.get(name, {}))
        absent = sorted(required - set(values))
        if absent:
            missing[name] = absent
        resolved[name] = values
    if missing:
        raise KeyError(f"The clustering execution contract is incomplete: {missing}.")
    return (
        resolved["neighbors_args"],
        resolved["umap_args"],
        resolved["leiden_args"],
    )


def _h5ad_safe_arguments(value: Any) -> Any:
    if value is None:
        return "null"
    if isinstance(value, dict):
        return {str(key): _h5ad_safe_arguments(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_h5ad_safe_arguments(item) for item in value]
    return value


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cl = cfg["clustering"]
    neighbors_args, umap_args, leiden_args = _required_clustering_arguments(cl)
    if "n_jobs" not in cl or int(cl["n_jobs"]) < 1:
        raise ValueError("clustering.n_jobs must be an explicit positive integer.")
    n_jobs = int(cl["n_jobs"])
    seed = int(cfg.get("run", {}).get("seed", 42))

    # PyNNDescent 0.5.x contains an internal n_jobs=-1 tree-building call.
    # Constrain joblib before importing Scanpy so the configured worker limit
    # also applies to that internal call on Windows.
    os.environ["LOKY_MAX_CPU_COUNT"] = str(n_jobs)
    import scanpy as sc

    sc.settings.n_jobs = n_jobs

    started_at = datetime.now(timezone.utc).isoformat()
    print(f"[cluster] Loading {args.inp}", flush=True)
    adata = sc.read_h5ad(args.inp)
    input_obs_sha256 = _index_sha256(adata.obs_names)
    input_var_sha256 = _index_sha256(adata.var_names)

    print(f"[cluster] Computing neighbors with {neighbors_args}", flush=True)
    sc.pp.neighbors(
        adata,
        use_rep=cl["use_rep"],
        random_state=seed,
        **neighbors_args,
    )
    print(f"[cluster] Computing UMAP with {umap_args}", flush=True)
    sc.tl.umap(adata, random_state=seed, **umap_args)
    print(f"[cluster] Computing Leiden clusters with {leiden_args}", flush=True)
    sc.tl.leiden(
        adata,
        resolution=float(cl["leiden_resolution"]),
        random_state=seed,
        **leiden_args,
    )
    adata.uns["clustering_provenance"] = {
        "started_at_utc": started_at,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(Path(args.inp).resolve()),
        "input_shape": [int(adata.n_obs), int(adata.n_vars)],
        "input_obs_names_sha256": input_obs_sha256,
        "input_var_names_sha256": input_var_sha256,
        "config_path": str(Path(args.config).resolve()),
        "resolved_config_sha256": _resolved_config_sha256(cfg),
        "seed": seed,
        "n_jobs": n_jobs,
        "use_rep": str(cl["use_rep"]),
        "leiden_resolution": float(cl["leiden_resolution"]),
        "neighbors_args": _h5ad_safe_arguments(neighbors_args),
        "umap_args": _h5ad_safe_arguments(umap_args),
        "leiden_args": _h5ad_safe_arguments(leiden_args),
        "scanpy_version": str(version("scanpy")),
    }

    print(f"[cluster] Writing {args.out}", flush=True)
    adata.write_h5ad(args.out)
    print("[cluster] Completed", flush=True)


if __name__ == "__main__":
    main()
