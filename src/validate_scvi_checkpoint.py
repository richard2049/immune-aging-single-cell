from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np

from .scvi_train import _index_sha256, _required_arguments, _resolved_config_sha256
from .utils import ensure_dir, load_config


def _sampled_file_fingerprint(path: Path, chunk_size: int = 1024 * 1024) -> dict[str, Any]:
    size = path.stat().st_size
    digest = hashlib.sha256()
    offsets = sorted({0, max(0, size // 2 - chunk_size // 2), max(0, size - chunk_size)})
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            digest.update(offset.to_bytes(8, "little"))
            digest.update(handle.read(chunk_size))
    return {
        "path": str(path.resolve()),
        "size_bytes": int(size),
        "sampled_sha256": digest.hexdigest(),
    }


def _record(
    checks: list[dict[str, Any]],
    name: str,
    passed: bool,
    detail: str,
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def _normalize_serialized(value: Any) -> Any:
    if isinstance(value, dict):
        normalized = {str(key): _normalize_serialized(item) for key, item in value.items()}
        return {key: item for key, item in normalized.items() if item is not None}
    if isinstance(value, np.ndarray):
        return [_normalize_serialized(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_normalize_serialized(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def validate(
    *,
    config_path: Path,
    input_path: Path,
    checkpoint_path: Path,
    report_path: Path,
) -> dict[str, Any]:
    cfg = load_config(config_path)
    scvicfg = cfg["scvi"]
    expected_model_args, expected_training_args = _required_arguments(scvicfg)
    model_dir = Path(str(cfg["project"]["out_dir"])) / "models" / "scvi_model"
    checks: list[dict[str, Any]] = []

    if not input_path.is_file() or not checkpoint_path.is_file():
        raise FileNotFoundError("Both the immutable input and scVI checkpoint are required.")

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

        latent = candidate.obsm.get("X_scVI")
        expected_latent = int(expected_model_args["n_latent"])
        latent_shape = tuple(latent.shape) if latent is not None else None
        _record(
            checks,
            "latent_shape",
            latent_shape == (candidate.n_obs, expected_latent),
            f"expected={(candidate.n_obs, expected_latent)}; observed={latent_shape}",
        )
        finite_latent = bool(latent is not None and np.isfinite(np.asarray(latent)).all())
        _record(checks, "latent_values_finite", finite_latent, f"observed={finite_latent}")

        provenance = candidate.uns.get("scvi_training_provenance")
        provenance_ok = isinstance(provenance, dict)
        _record(checks, "provenance_present", provenance_ok, f"type={type(provenance).__name__}")
        if provenance_ok:
            expected = {
                "seed": int(cfg.get("run", {}).get("seed", 42)),
                "layer": "X"
                if scvicfg.get("layer") in (None, "", "null", "None")
                else str(scvicfg["layer"]),
                "categorical_covariates": list(scvicfg.get("categorical_covariates", [])),
                "continuous_covariates": list(scvicfg.get("continuous_covariates", [])),
                "allow_missing_covariates": bool(scvicfg.get("allow_missing_covariates", False)),
                "model_args": expected_model_args,
                "training_args": expected_training_args,
                "accelerator": str(scvicfg["accelerator"]),
                "devices": int(scvicfg.get("devices", 1)),
                "input_shape": [int(source.n_obs), int(source.n_vars)],
                "input_obs_names_sha256": input_obs_hash,
                "input_var_names_sha256": input_var_hash,
                "resolved_config_sha256": _resolved_config_sha256(cfg),
            }
            for key, value in expected.items():
                observed = provenance.get(key)
                passed = _normalize_serialized(observed) == _normalize_serialized(value)
                _record(
                    checks,
                    f"provenance_matches__{key}",
                    passed,
                    f"expected={value!r}; observed={observed!r}",
                )
            for key in (
                "started_at_utc",
                "completed_at_utc",
                "scvi_tools_version",
                "torch_version",
                "anndata_version",
                "scanpy_version",
            ):
                value = str(provenance.get(key, "")).strip()
                _record(checks, f"provenance_nonempty__{key}", bool(value), value or "missing")

        model_files = (
            sorted(path for path in model_dir.rglob("*") if path.is_file())
            if model_dir.is_dir()
            else []
        )
        model_bytes = sum(path.stat().st_size for path in model_files)
        _record(
            checks,
            "model_persisted",
            bool(model_files) and model_bytes > 0,
            f"directory={model_dir}; files={len(model_files)}; bytes={model_bytes}",
        )
    finally:
        source.file.close()
        candidate.file.close()

    failed = [check["name"] for check in checks if not check["passed"]]
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if not failed else "failed",
        "failed_checks": failed,
        "checks": checks,
        "input": _sampled_file_fingerprint(input_path),
        "checkpoint": _sampled_file_fingerprint(checkpoint_path),
        "model_directory": str(model_dir.resolve()),
        "config_path": str(config_path.resolve()),
    }
    ensure_dir(report_path.parent)
    temporary = report_path.with_name(report_path.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2), encoding="utf-8")
    os.replace(temporary, report_path)
    if failed:
        raise ValueError("scVI checkpoint validation failed: " + ", ".join(failed))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a reconstructed scVI checkpoint.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    validate(
        config_path=Path(args.config),
        input_path=Path(args.inp),
        checkpoint_path=Path(args.checkpoint),
        report_path=Path(args.report),
    )


if __name__ == "__main__":
    main()
