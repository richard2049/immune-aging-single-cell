from __future__ import annotations
import argparse
import os
from pathlib import Path
import shutil

import anndata as ad
import numpy as np
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
        adata.write_h5ad(args.out)
    elif dataset == "custom_h5ad":
        in_path = Path(cfg["paths"]["input_h5ad"]).resolve()
        out_path = Path(args.out).resolve()

        if not in_path.exists():
            raise FileNotFoundError(f"Input h5ad not found: {in_path}")

        # Optional pilot mode for very large datasets.
        # If run.max_cells is set, sample cells from backed AnnData without
        # loading the full matrix into memory.
        max_cells = cfg.get("run", {}).get("max_cells")
        if max_cells is not None:
            max_cells = int(max_cells)
            if max_cells <= 0:
                raise ValueError("run.max_cells must be > 0 when provided")

            seed = int(cfg.get("run", {}).get("seed", 0))
            drop_raw_on_subset = bool(cfg.get("run", {}).get("drop_raw_on_subset", True))
            adata_backed = ad.read_h5ad(str(in_path), backed="r")
            try:
                # Backed `.to_memory()` also materializes `.raw` when present.
                # For large files this can dominate RAM and trigger OOM during
                # subsampling; keep it optional but enabled by default.
                if drop_raw_on_subset and adata_backed.raw is not None:
                    adata_backed.raw = None

                n_obs = int(adata_backed.n_obs)

                if n_obs > max_cells:
                    rng = np.random.default_rng(seed)
                    idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))
                    adata = adata_backed[idx, :].to_memory()
                else:
                    adata = adata_backed.to_memory()
            finally:
                if hasattr(adata_backed, "file") and getattr(adata_backed, "file", None) is not None:
                    adata_backed.file.close()

            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            adata.write_h5ad(args.out)
            return

        # Default for large real datasets: avoid in-memory reload/rewrite.
        # Hard-link when possible; fall back to file copy.
        ensure_dir(out_path.parent)
        if out_path.exists():
            out_path.unlink()

        try:
            os.link(in_path, out_path)
        except OSError:
            shutil.copy2(in_path, out_path)
        return
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

if __name__ == "__main__":
    main()
