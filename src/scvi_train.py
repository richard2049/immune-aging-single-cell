from __future__ import annotations

import argparse
import scanpy as sc
import scipy.sparse as sp
import scvi

from .utils import load_config, ensure_dir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    scvicfg = cfg["scvi"]
    out_dir = cfg["project"]["out_dir"]

    ensure_dir(f"{out_dir}/models")

    adata = sc.read_h5ad(args.inp)

    # Ensure sparse counts layer exists
    if "counts" not in adata.layers:
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(adata.X)
        adata.layers["counts"] = adata.X.copy()

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=scvicfg["layer"],
        categorical_covariate_keys=scvicfg["categorical_covariates"],
        continuous_covariate_keys=scvicfg["continuous_covariates"],
    )

    model = scvi.model.SCVI(adata)

    if scvicfg["accelerator"] == "gpu":
        model.train(accelerator="gpu", devices=int(scvicfg["devices"]))
    else:
        model.train(accelerator="cpu")

    adata.obsm["X_scVI"] = model.get_latent_representation()

    # Save model weights (useful for reproducibility)
    model.save(f"{out_dir}/models/scvi_model", overwrite=True)

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()