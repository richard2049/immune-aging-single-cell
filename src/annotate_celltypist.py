from __future__ import annotations

import argparse
import scanpy as sc
import celltypist
import numpy as np

from .utils import load_config


def _to_1d_series(x):
    """Convert pandas Series/DataFrame to a 1D Series."""
    ndim = getattr(x, "ndim", None)
    if ndim == 1:
        return x
    if ndim == 2:
        cols = list(getattr(x, "columns", []))
        # Common column names across versions
        for c in ("majority_voting", "predicted_labels", "label", "cell_type"):
            if c in cols:
                return x[c]
        return x.iloc[:, 0]
    return x.squeeze()


def _compute_confidence(res):
    """
    Try to compute a per-cell confidence score.
    - Prefer probability_matrix (max prob per cell)
    - Fallback to decision_matrix (max score per cell)
    - Otherwise return None
    """
    pm = getattr(res, "probability_matrix", None)
    if pm is not None:
        try:
            return pm.max(axis=1)
        except (TypeError, ValueError, AttributeError):
            pass

    dm = getattr(res, "decision_matrix", None)
    if dm is not None:
        try:
            return dm.max(axis=1)
        except (TypeError, ValueError, AttributeError):
            pass

    return None


def _extract_labels_and_confidence(res, obs_names):
    labels_raw = getattr(res, "majority_voting", None)
    if labels_raw is None:
        labels_raw = res.predicted_labels

    labels = _to_1d_series(labels_raw)
    if hasattr(labels, "reindex"):
        labels = labels.reindex(obs_names)
    labels = labels.astype(str)

    conf = _compute_confidence(res)
    if conf is not None:
        conf = _to_1d_series(conf)
        if hasattr(conf, "reindex"):
            conf = conf.reindex(obs_names)
        try:
            conf = conf.astype(float)
        except (TypeError, ValueError):
            conf = None

    return labels, conf


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ctcfg = cfg["celltypist"]

    adata = sc.read_h5ad(args.inp)

    if not bool(ctcfg.get("enabled", True)):
        adata.write_h5ad(args.out)
        return

    chunk_size = int(ctcfg.get("chunk_size", 0))
    use_chunking = chunk_size > 0 and adata.n_obs > chunk_size
    if use_chunking and bool(ctcfg.get("majority_voting", True)):
        print(
            "CellTypist chunk mode enabled; overriding majority_voting=False "
            "to avoid chunk-wise voting artifacts."
        )

    majority_voting = bool(ctcfg.get("majority_voting", True)) and not use_chunking

    # CellTypist expects log-normalized expression
    target_sum = float(ctcfg.get("normalize_target_sum", 1e4))

    if use_chunking:
        labels_out = np.empty(adata.n_obs, dtype=object)
        conf_out = np.full(adata.n_obs, np.nan, dtype=float)

        for start in range(0, adata.n_obs, chunk_size):
            end = min(start + chunk_size, adata.n_obs)
            print(f"Annotating cells {start}:{end} / {adata.n_obs}")

            tmp = adata[start:end].copy()
            sc.pp.normalize_total(tmp, target_sum=target_sum)
            sc.pp.log1p(tmp)

            res = celltypist.annotate(
                tmp,
                model=ctcfg["model"],
                majority_voting=majority_voting,
            )

            labels, conf = _extract_labels_and_confidence(res, tmp.obs_names)
            labels_out[start:end] = labels.to_numpy()
            if conf is not None:
                conf_out[start:end] = conf.to_numpy()

        adata.obs["cell_type"] = labels_out
        if np.isfinite(conf_out).any():
            adata.obs["cell_type_confidence"] = conf_out
    else:
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=target_sum)
        sc.pp.log1p(tmp)

        res = celltypist.annotate(
            tmp,
            model=ctcfg["model"],
            majority_voting=majority_voting,
        )

        labels, conf = _extract_labels_and_confidence(res, adata.obs_names)
        adata.obs["cell_type"] = labels.to_numpy()
        if conf is not None:
            adata.obs["cell_type_confidence"] = conf.to_numpy()

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
