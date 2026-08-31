from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from .cluster import _h5ad_safe_arguments
from .scvi_train import _index_sha256
from .utils import ensure_dir, load_config
from .validate_scvi_checkpoint import _normalize_serialized, _sampled_file_fingerprint


def _record(
    checks: list[dict[str, Any]],
    name: str,
    passed: bool,
    detail: str,
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def _expected_provenance(config: dict[str, Any]) -> dict[str, Any]:
    clustering = config["clustering"]
    return {
        "seed": int(config.get("run", {}).get("seed", 42)),
        "n_jobs": int(clustering["n_jobs"]),
        "use_rep": str(clustering["use_rep"]),
        "leiden_resolution": float(clustering["leiden_resolution"]),
        "neighbors_args": _h5ad_safe_arguments(clustering["neighbors_args"]),
        "umap_args": _h5ad_safe_arguments(clustering["umap_args"]),
        "leiden_args": _h5ad_safe_arguments(clustering["leiden_args"]),
    }


def validate(
    *,
    config_path: Path,
    input_path: Path,
    checkpoint_path: Path,
    report_path: Path,
    historical_path: Path | None = None,
) -> dict[str, Any]:
    if not input_path.is_file() or not checkpoint_path.is_file():
        raise FileNotFoundError(
            "Both the validated scVI input and clustered checkpoint are required."
        )

    config = load_config(config_path)
    expected_provenance = _expected_provenance(config)
    expected_umap_components = int(config["clustering"]["umap_args"]["n_components"])
    checks: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}

    source = ad.read_h5ad(input_path, backed="r")
    candidate = ad.read_h5ad(checkpoint_path, backed="r")
    try:
        _record(
            checks,
            "shape_preserved",
            source.shape == candidate.shape,
            f"input={source.shape}; checkpoint={candidate.shape}",
        )
        _record(
            checks,
            "cell_identifiers_unique",
            bool(candidate.obs_names.is_unique),
            f"n_cells={candidate.n_obs}",
        )
        _record(
            checks,
            "gene_identifiers_unique",
            bool(candidate.var_names.is_unique),
            f"n_genes={candidate.n_vars}",
        )
        input_obs_hash = _index_sha256(source.obs_names)
        input_var_hash = _index_sha256(source.var_names)
        output_obs_hash = _index_sha256(candidate.obs_names)
        output_var_hash = _index_sha256(candidate.var_names)
        _record(
            checks,
            "cell_identity_and_order_preserved",
            input_obs_hash == output_obs_hash,
            f"input={input_obs_hash}; checkpoint={output_obs_hash}",
        )
        _record(
            checks,
            "gene_identity_and_order_preserved",
            input_var_hash == output_var_hash,
            f"input={input_var_hash}; checkpoint={output_var_hash}",
        )

        umap = candidate.obsm.get("X_umap")
        umap_shape = tuple(umap.shape) if umap is not None else None
        expected_umap_shape = (candidate.n_obs, expected_umap_components)
        _record(
            checks,
            "umap_shape",
            umap_shape == expected_umap_shape,
            f"expected={expected_umap_shape}; observed={umap_shape}",
        )
        umap_finite = bool(umap is not None and np.isfinite(np.asarray(umap)).all())
        _record(checks, "umap_values_finite", umap_finite, f"observed={umap_finite}")

        for key in ("distances", "connectivities"):
            graph = candidate.obsp.get(key)
            graph_shape = tuple(graph.shape) if graph is not None else None
            expected_graph_shape = (candidate.n_obs, candidate.n_obs)
            graph_nnz = int(graph.nnz) if graph is not None else 0
            graph_finite = bool(
                graph is not None and graph_nnz > 0 and np.isfinite(graph.data).all()
            )
            _record(
                checks,
                f"{key}_graph",
                graph_shape == expected_graph_shape and graph_finite,
                f"shape={graph_shape}; nnz={graph_nnz}; finite={graph_finite}",
            )

        leiden = candidate.obs.get("leiden")
        leiden_missing = int(leiden.isna().sum()) if leiden is not None else candidate.n_obs
        cluster_sizes = leiden.astype(str).value_counts() if leiden is not None else None
        n_clusters = int(len(cluster_sizes)) if cluster_sizes is not None else 0
        _record(
            checks,
            "leiden_labels_complete",
            leiden is not None and leiden_missing == 0 and n_clusters > 0,
            f"missing={leiden_missing}; n_clusters={n_clusters}",
        )
        summary["candidate_clusters"] = n_clusters
        summary["candidate_cluster_sizes"] = (
            {str(key): int(value) for key, value in cluster_sizes.items()}
            if cluster_sizes is not None
            else {}
        )

        provenance = candidate.uns.get("clustering_provenance")
        _record(
            checks,
            "provenance_present",
            isinstance(provenance, dict),
            f"type={type(provenance).__name__}",
        )
        if isinstance(provenance, dict):
            for key, expected in expected_provenance.items():
                observed = provenance.get(key)
                passed = _normalize_serialized(observed) == _normalize_serialized(expected)
                _record(
                    checks,
                    f"provenance_matches__{key}",
                    passed,
                    f"expected={expected!r}; observed={observed!r}",
                )
            for key in (
                "started_at_utc",
                "completed_at_utc",
                "input_path",
                "input_obs_names_sha256",
                "input_var_names_sha256",
                "resolved_config_sha256",
                "scanpy_version",
            ):
                value = str(provenance.get(key, "")).strip()
                _record(checks, f"provenance_nonempty__{key}", bool(value), value or "missing")
            _record(
                checks,
                "provenance_input_identity",
                provenance.get("input_obs_names_sha256") == input_obs_hash
                and provenance.get("input_var_names_sha256") == input_var_hash,
                "ordered input hashes compared",
            )
    finally:
        source.file.close()
        candidate.file.close()

    if historical_path is not None and historical_path.is_file():
        historical = ad.read_h5ad(historical_path, backed="r")
        candidate = ad.read_h5ad(checkpoint_path, backed="r")
        try:
            historical_hash = _index_sha256(historical.obs_names)
            candidate_hash = _index_sha256(candidate.obs_names)
            if (
                historical_hash == candidate_hash
                and "leiden" in historical.obs
                and "leiden" in candidate.obs
            ):
                historical_labels = historical.obs["leiden"].astype(str)
                candidate_labels = candidate.obs["leiden"].astype(str)
                summary["historical_comparison"] = {
                    "historical_clusters": int(historical_labels.nunique()),
                    "adjusted_rand_index": float(
                        adjusted_rand_score(historical_labels, candidate_labels)
                    ),
                    "normalized_mutual_information": float(
                        normalized_mutual_info_score(historical_labels, candidate_labels)
                    ),
                    "acceptance_role": "informative_only",
                }
            else:
                summary["historical_comparison"] = {
                    "available": False,
                    "reason": "cell identity/order or Leiden labels are incompatible",
                }
        finally:
            historical.file.close()
            candidate.file.close()

    failed = [check["name"] for check in checks if not check["passed"]]
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if not failed else "failed",
        "failed_checks": failed,
        "checks": checks,
        "summary": summary,
        "input": _sampled_file_fingerprint(input_path),
        "checkpoint": _sampled_file_fingerprint(checkpoint_path),
        "config_path": str(config_path.resolve()),
        "historical_path": str(historical_path.resolve()) if historical_path else None,
    }
    ensure_dir(report_path.parent)
    temporary = report_path.with_name(report_path.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2), encoding="utf-8")
    os.replace(temporary, report_path)
    if failed:
        raise ValueError("Clustered checkpoint validation failed: " + ", ".join(failed))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a clustered AnnData checkpoint.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--historical")
    args = parser.parse_args()
    validate(
        config_path=Path(args.config),
        input_path=Path(args.inp),
        checkpoint_path=Path(args.checkpoint),
        report_path=Path(args.report),
        historical_path=Path(args.historical) if args.historical else None,
    )


if __name__ == "__main__":
    main()
