from __future__ import annotations

import argparse
import scanpy as sc
import celltypist

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
        except Exception:
            pass

    dm = getattr(res, "decision_matrix", None)
    if dm is not None:
        try:
            return dm.max(axis=1)
        except Exception:
            pass

    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ctcfg = cfg["celltypist"]

    adata = sc.read_h5ad(args.inp)

    # CellTypist expects log-normalized expression
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

    res = celltypist.annotate(
        tmp,
        model=ctcfg["model"],
        majority_voting=bool(ctcfg.get("majority_voting", True)),
    )

    # Labels: prefer majority voting labels if available
    labels_raw = getattr(res, "majority_voting", None)
    if labels_raw is None:
        labels_raw = res.predicted_labels

    labels = _to_1d_series(labels_raw)
    if hasattr(labels, "reindex"):
        labels = labels.reindex(adata.obs_names)

    adata.obs["cell_type"] = labels.astype(str).to_numpy()

    conf = _compute_confidence(res)
    if conf is not None:
        conf = _to_1d_series(conf)
        if hasattr(conf, "reindex"):
            conf = conf.reindex(adata.obs_names)
        # store as float
        try:
            adata.obs["cell_type_confidence"] = conf.astype(float).to_numpy()
        except Exception:
            pass

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
