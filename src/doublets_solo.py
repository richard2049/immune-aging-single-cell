from __future__ import annotations

import argparse
import scanpy as sc
import scipy.sparse as sp
import scvi

from .utils import load_config


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    dbl = cfg["doublets"]
    scvicfg = cfg["scvi"]

    adata = sc.read_h5ad(args.inp)

    # If disabled, pass-through
    if not bool(dbl.get("enabled", True)):
        adata.write_h5ad(args.out)
        return

    # Guardrail: SOLO is unstable on tiny datasets
    if adata.n_obs < 200:
        adata.write_h5ad(args.out)
        return

    # Train SOLO on a filtered HVG subset (lighter & common practice)
    train = adata.copy()
    if not sp.issparse(train.X):
        train.X = sp.csr_matrix(train.X)

    # Prevent min_cells from removing all genes when cell count is small
    min_cells_cfg = int(dbl.get("filter_genes_min_cells", 10))
    min_cells = min(min_cells_cfg, max(1, train.n_obs - 1))

    sc.pp.filter_genes(train, min_cells=min_cells)
    if train.n_vars == 0 or train.n_obs == 0:
        adata.write_h5ad(args.out)
        return

    hvg_flavor = str(dbl.get("hvg_flavor", "seurat_v3"))

    sc.pp.highly_variable_genes(
        train,
        n_top_genes=int(dbl.get("hvg_top", 2000)),
        subset=True,
        flavor=hvg_flavor,
    )

    if train.n_vars == 0:
        adata.write_h5ad(args.out)
        return

    scvi.model.SCVI.setup_anndata(train)
    vae = scvi.model.SCVI(train)

    if str(scvicfg.get("accelerator", "cpu")) == "gpu":
        vae.train(accelerator="gpu", devices=int(scvicfg.get("devices", 1)))
        solo = scvi.external.SOLO.from_scvi_model(vae)
        solo.train(accelerator="gpu", devices=int(scvicfg.get("devices", 1)))
    else:
        vae.train(accelerator="cpu")
        solo = scvi.external.SOLO.from_scvi_model(vae)
        solo.train(accelerator="cpu")

    df = solo.predict()
    df["prediction"] = solo.predict(soft=False)
    df["dif"] = df["doublet"] - df["singlet"]

    dif_thr = float(dbl.get("dif_threshold", 0.5))
    doublets = df[(df["prediction"] == "doublet") & (df["dif"] > dif_thr)]

    adata.obs["doublet"] = adata.obs_names.isin(doublets.index)
    adata = adata[~adata.obs["doublet"]].copy()

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()