from __future__ import annotations

import argparse
import scanpy as sc

from .utils import load_config


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cl = cfg["clustering"]

    adata = sc.read_h5ad(args.inp)

    sc.pp.neighbors(adata, use_rep=cl["use_rep"])
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=float(cl["leiden_resolution"]))

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()