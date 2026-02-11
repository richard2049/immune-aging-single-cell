from __future__ import annotations

import argparse
import numpy as np
import scanpy as sc
import scipy.sparse as sp

from .utils import load_config


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    qc = cfg["qc"]

    adata = sc.read_h5ad(args.inp)

    # Ensure sparse counts + keep raw counts in layers["counts"]
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.layers["counts"] = adata.X.copy()

    # mt + ribo flags (portable defaults)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPL", "RPS"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Filter genes/cells
    sc.pp.filter_genes(adata, min_cells=int(qc["min_cells_per_gene"]))
    sc.pp.filter_cells(adata, min_genes=int(qc["min_genes_per_cell"]))

    upper_lim = np.quantile(
        adata.obs["n_genes_by_counts"].to_numpy(),
        float(qc["n_genes_upper_quantile"]),
    )

    adata = adata[adata.obs["n_genes_by_counts"] < upper_lim].copy()
    adata = adata[adata.obs["pct_counts_mt"] < float(qc["max_pct_counts_mt"])].copy()
    adata = adata[adata.obs["pct_counts_ribo"] < float(qc["max_pct_counts_ribo"])].copy()

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()