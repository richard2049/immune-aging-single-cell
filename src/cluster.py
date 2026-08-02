from __future__ import annotations

import argparse

from .utils import load_config


def main() -> None:
    import scanpy as sc

    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cl = cfg["clustering"]
    seed = int(cfg.get("run", {}).get("seed", 42))

    adata = sc.read_h5ad(args.inp)

    sc.pp.neighbors(adata, use_rep=cl["use_rep"], random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=float(cl["leiden_resolution"]),
        random_state=seed,
    )
    adata.uns["clustering_provenance"] = {
        "seed": seed,
        "use_rep": str(cl["use_rep"]),
        "leiden_resolution": float(cl["leiden_resolution"]),
    }

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
