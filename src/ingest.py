from __future__ import annotations
import argparse
import scanpy as sc
from .utils import load_config, ensure_dir

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ensure_dir(cfg["project"]["out_dir"])

    dataset = cfg["run"]["dataset"]

    if dataset == "demo_pbmc3k":
        adata = sc.datasets.pbmc3k()
        # ensure gene symbols
        adata.var_names_make_unique()
    elif dataset == "custom_h5ad":
        adata = sc.read_h5ad(cfg["paths"]["input_h5ad"])
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    adata.write_h5ad(args.out)

if __name__ == "__main__":
    main()