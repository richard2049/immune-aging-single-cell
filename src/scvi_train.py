from __future__ import annotations

import argparse

from .utils import load_config, ensure_dir


def _resolve_covariates(obs_columns, scvicfg: dict) -> tuple[list[str], list[str]]:
    raw_cat = [str(key) for key in scvicfg.get("categorical_covariates", [])]
    raw_cont = [str(key) for key in scvicfg.get("continuous_covariates", [])]
    cat_keys = [key for key in raw_cat if key in obs_columns]
    cont_keys = [key for key in raw_cont if key in obs_columns]
    missing = {
        "categorical": sorted(set(raw_cat) - set(cat_keys)),
        "continuous": sorted(set(raw_cont) - set(cont_keys)),
    }
    if any(missing.values()) and not bool(
        scvicfg.get("allow_missing_covariates", False)
    ):
        raise KeyError(
            "Configured scVI covariates are missing from adata.obs: "
            f"{missing}. Set scvi.allow_missing_covariates=true only for an "
            "explicitly documented alternative model."
        )
    if any(missing.values()):
        print(f"[scvi_train] Explicitly allowing missing covariates: {missing}")
    return cat_keys, cont_keys


def main() -> None:
    import scanpy as sc
    import scipy.sparse as sp
    import scvi

    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    scvicfg = cfg["scvi"]
    out_dir = cfg["project"]["out_dir"]
    seed = int(cfg.get("run", {}).get("seed", 42))
    scvi.settings.seed = seed
    layer = scvicfg.get("layer")
    if layer in ("", "null", "None"):
        layer = None

    ensure_dir(f"{out_dir}/models")

    adata = sc.read_h5ad(args.inp)

    # Ensure sparse representation for memory efficiency.
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    # Create a dedicated layer only when explicitly requested.
    if layer is not None and layer not in adata.layers:
        adata.layers[layer] = adata.X.copy()

    cat_keys, cont_keys = _resolve_covariates(adata.obs.columns, scvicfg)

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        categorical_covariate_keys=cat_keys,
        continuous_covariate_keys=cont_keys,
    )

    model = scvi.model.SCVI(adata)

    if scvicfg["accelerator"] == "gpu":
        model.train(accelerator="gpu", devices=int(scvicfg["devices"]))
    else:
        model.train(accelerator="cpu")

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.uns["scvi_training_provenance"] = {
        "seed": seed,
        "layer": "X" if layer is None else str(layer),
        "categorical_covariates": cat_keys,
        "continuous_covariates": cont_keys,
        "allow_missing_covariates": bool(
            scvicfg.get("allow_missing_covariates", False)
        ),
    }

    # Save model weights (useful for reproducibility)
    model.save(f"{out_dir}/models/scvi_model", overwrite=True)

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
