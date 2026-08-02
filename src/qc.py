from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp

from .utils import ensure_dir, load_config


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    qc = cfg["qc"]

    # Large-data passthrough mode: skip in-memory QC computations.
    if not bool(qc.get("enabled", True)):
        in_path = Path(args.inp).resolve()
        out_path = Path(args.out).resolve()
        ensure_dir(out_path.parent)
        if out_path.exists():
            out_path.unlink()
        try:
            os.link(in_path, out_path)
        except OSError:
            shutil.copy2(in_path, out_path)
        return

    adata = sc.read_h5ad(args.inp)

    # Ensure sparse matrix representation before QC operations.
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

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

    obs_mask = (
        (adata.obs["n_genes_by_counts"] < upper_lim)
        & (adata.obs["pct_counts_mt"] < float(qc["max_pct_counts_mt"]))
        & (adata.obs["pct_counts_ribo"] < float(qc["max_pct_counts_ribo"]))
    )
    adata = adata[obs_mask].copy()

    # Store filtered raw counts for downstream scVI setup.
    adata.layers["counts"] = adata.X.copy()

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
