from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

from .utils import ensure_dir, load_config


def _resolve_covariates(
    obs_columns: Any,
    scvicfg: dict[str, Any],
) -> tuple[list[str], list[str]]:
    raw_cat = [str(key) for key in scvicfg.get("categorical_covariates", [])]
    raw_cont = [str(key) for key in scvicfg.get("continuous_covariates", [])]
    cat_keys = [key for key in raw_cat if key in obs_columns]
    cont_keys = [key for key in raw_cont if key in obs_columns]
    missing = {
        "categorical": sorted(set(raw_cat) - set(cat_keys)),
        "continuous": sorted(set(raw_cont) - set(cont_keys)),
    }
    if any(missing.values()) and not bool(scvicfg.get("allow_missing_covariates", False)):
        raise KeyError(
            "Configured scVI covariates are missing from adata.obs: "
            f"{missing}. Set scvi.allow_missing_covariates=true only for an "
            "explicitly documented alternative model."
        )
    if any(missing.values()):
        print(f"[scvi_train] Explicitly allowing missing covariates: {missing}")
    return cat_keys, cont_keys


def _required_arguments(
    scvicfg: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    model_args = dict(scvicfg.get("model_args", {}))
    training_args = dict(scvicfg.get("training_args", {}))
    required_model = {
        "n_hidden",
        "n_latent",
        "n_layers",
        "dropout_rate",
        "dispersion",
        "gene_likelihood",
        "use_observed_lib_size",
        "latent_distribution",
    }
    required_training = {
        "max_epochs",
        "batch_size",
        "early_stopping",
        "train_size",
        "validation_size",
        "shuffle_set_split",
        "load_sparse_tensor",
    }
    missing_model = sorted(required_model - set(model_args))
    missing_training = sorted(required_training - set(training_args))
    if missing_model or missing_training:
        raise KeyError(
            "The scVI execution contract must explicitly define all model and "
            "training arguments. Missing model_args="
            f"{missing_model}; training_args={missing_training}."
        )
    return model_args, training_args


def _index_sha256(values: Any) -> str:
    digest = hashlib.sha256()
    for value in values:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
    return digest.hexdigest()


def _resolved_config_sha256(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _software_versions(torch_version: Any) -> dict[str, str]:
    """Return plain strings that AnnData can serialize without type loss."""
    return {
        "scvi_tools_version": str(version("scvi-tools")),
        "torch_version": str(torch_version),
        "anndata_version": str(version("anndata")),
        "scanpy_version": str(version("scanpy")),
    }


def main() -> None:
    import scanpy as sc
    import scipy.sparse as sp
    import scvi
    import torch

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    scvicfg = cfg["scvi"]
    model_args, training_args = _required_arguments(scvicfg)
    out_dir = Path(str(cfg["project"]["out_dir"]))
    model_dir = out_dir / "models" / "scvi_model"
    seed = int(cfg.get("run", {}).get("seed", 42))
    scvi.settings.seed = seed
    layer = scvicfg.get("layer")
    if layer in ("", "null", "None"):
        layer = None

    ensure_dir(out_dir / "models")
    started_at = datetime.now(timezone.utc).isoformat()
    print(f"[scvi_train] Loading {Path(args.inp).resolve()}", flush=True)
    adata = sc.read_h5ad(args.inp)
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("scVI input requires unique cell and gene identifiers.")
    input_obs_sha256 = _index_sha256(adata.obs_names)
    input_var_sha256 = _index_sha256(adata.var_names)

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    if layer is not None and layer not in adata.layers:
        adata.layers[layer] = adata.X.copy()

    cat_keys, cont_keys = _resolve_covariates(adata.obs.columns, scvicfg)
    print(
        f"[scvi_train] Setting up AnnData: cells={adata.n_obs}, genes={adata.n_vars}, "
        f"categorical_covariates={cat_keys}, continuous_covariates={cont_keys}",
        flush=True,
    )
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        categorical_covariate_keys=cat_keys,
        continuous_covariate_keys=cont_keys,
    )

    print(f"[scvi_train] Initializing model with {model_args}", flush=True)
    model = scvi.model.SCVI(adata, **model_args)
    accelerator = str(scvicfg["accelerator"])
    devices = int(scvicfg.get("devices", 1))
    print(
        f"[scvi_train] Training on accelerator={accelerator}, devices={devices}, "
        f"training_args={training_args}",
        flush=True,
    )
    model.train(
        accelerator=accelerator,
        devices=devices,
        **training_args,
    )

    print("[scvi_train] Extracting latent representation", flush=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.uns["scvi_training_provenance"] = {
        "started_at_utc": started_at,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(Path(args.inp).resolve()),
        "input_shape": [int(adata.n_obs), int(adata.n_vars)],
        "input_obs_names_sha256": input_obs_sha256,
        "input_var_names_sha256": input_var_sha256,
        "config_path": str(Path(args.config).resolve()),
        "resolved_config_sha256": _resolved_config_sha256(cfg),
        "seed": seed,
        "layer": "X" if layer is None else str(layer),
        "categorical_covariates": cat_keys,
        "continuous_covariates": cont_keys,
        "allow_missing_covariates": bool(scvicfg.get("allow_missing_covariates", False)),
        "model_args": model_args,
        "training_args": training_args,
        "accelerator": accelerator,
        "devices": devices,
        **_software_versions(torch.__version__),
    }

    print(f"[scvi_train] Saving model to {model_dir.resolve()}", flush=True)
    model.save(model_dir, overwrite=True)
    print(f"[scvi_train] Writing checkpoint to {Path(args.out).resolve()}", flush=True)
    adata.write_h5ad(args.out)
    print("[scvi_train] Completed", flush=True)


if __name__ == "__main__":
    main()
