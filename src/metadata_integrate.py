from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import re

import pandas as pd
import scanpy as sc

from .utils import ensure_dir, load_config

STANDARD_ALIASES = {
    "cell_id": {"cellid", "cell", "cellbarcode", "barcode", "unnamed0"},
    "donor_id": {
        "donorid",
        "donor",
        "subjectid",
        "subject",
        "participantid",
        "participant",
        "patientid",
        "patient",
    },
    "age": {"age", "donorage", "ageyears", "ageyrs", "years"},
    "sex": {"sex", "gender"},
    "age_group": {"agegroup", "agegrp"},
    "batch": {"batch", "batchid"},
    "tube_id": {"tubeid"},
    "file_name": {"filename", "filepath"},
}


def _canon(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def _resolve_col(name: str, columns: pd.Index) -> str | None:
    if name in columns:
        return name
    target = _canon(name)
    for c in columns:
        if _canon(c) == target:
            return c
    return None


def _normalize_metadata_columns(meta: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rename_map = {}

    # User overrides first (exact source column names)
    user_map = cfg.get("rename_columns", {}) or {}
    for src, dst in user_map.items():
        if src in meta.columns and dst not in meta.columns:
            rename_map[src] = dst

    # Automatic normalization for common clinical metadata headers
    if bool(cfg.get("auto_standardize_columns", True)):
        norm_to_col = {}
        for c in meta.columns:
            n = _canon(c)
            if n not in norm_to_col:
                norm_to_col[n] = c

        existing_targets = set(meta.columns) | set(rename_map.values())
        for target, aliases in STANDARD_ALIASES.items():
            if target in existing_targets:
                continue
            source = None
            for alias in aliases:
                source = norm_to_col.get(alias)
                if source is not None:
                    break
            if source is not None and source not in rename_map:
                rename_map[source] = target

    if rename_map:
        meta = meta.rename(columns=rename_map)
    return meta


def _passthrough(inp: str, out: str) -> None:
    in_path = Path(inp).resolve()
    out_path = Path(out).resolve()
    ensure_dir(out_path.parent)
    if out_path.exists():
        out_path.unlink()
    try:
        os.link(in_path, out_path)
    except OSError:
        shutil.copy2(in_path, out_path)


def _add_ids_from_obs_names(adata, cfg: dict) -> None:
    """Parse sample / donor IDs from obs_names of the form SAMPLE_BARCODE."""
    field_sample = cfg.get("sample_field", "sample_id")
    field_barcode = cfg.get("barcode_field", "cell_barcode")
    field_donor = cfg.get("donor_field", "donor_id")
    field_batch = cfg.get("batch_field", "batch")
    sample_delim = str(cfg.get("sample_delim", "auto"))
    donor_index = int(cfg.get("donor_index", 1))

    obs_names = pd.Index(adata.obs_names.astype(str))
    sample_token = obs_names.to_series(index=obs_names).str.split("_", n=1).str[0]
    cell_barcode = obs_names.to_series(index=obs_names).str.split("_", n=1).str[-1]

    if sample_delim.lower() == "auto":
        candidates = ["-", ".", "|", ":"]
        best_delim = "-"
        best_hits = -1
        probe = sample_token.head(min(100000, sample_token.shape[0]))
        for d in candidates:
            hits = int(probe.str.contains(re.escape(d), regex=True).sum())
            if hits > best_hits:
                best_hits = hits
                best_delim = d
        sample_delim = best_delim

    adata.obs[field_sample] = sample_token.to_numpy()
    adata.obs[field_barcode] = cell_barcode.to_numpy()
    adata.obs[field_batch] = sample_token.to_numpy()

    sample_parts = sample_token.str.split(sample_delim, regex=False)
    donor = sample_parts.str[donor_index]
    adata.obs[field_donor] = donor.to_numpy()


def _merge_external_metadata(adata, cfg: dict) -> None:
    table_path = cfg.get("table_path")
    if not table_path:
        return

    path = Path(table_path)
    if not path.exists():
        raise FileNotFoundError(f"metadata.table_path not found: {path}")

    sep = str(cfg.get("table_sep", ","))
    obs_key = str(cfg.get("obs_join_key", "donor_id"))
    requested_table_key = str(cfg.get("table_join_key", obs_key))

    meta = pd.read_csv(path, sep=sep)
    meta = _normalize_metadata_columns(meta, cfg)

    table_key = _resolve_col(requested_table_key, meta.columns)
    if table_key is None:
        raise KeyError(f"metadata.table_join_key '{requested_table_key}' missing in table")

    use_obs_names = obs_key in {"obs_names", "obs_name", "__obs_names__", "cell_id"}
    if not use_obs_names and obs_key not in adata.obs.columns:
        raise KeyError(f"metadata.obs_join_key '{obs_key}' missing in adata.obs")

    cols = cfg.get("table_columns")
    if cols:
        keep = [table_key]
        for c in cols:
            resolved = _resolve_col(str(c), meta.columns)
            if resolved is not None and resolved != table_key and resolved not in keep:
                keep.append(resolved)
        meta = meta.loc[:, keep]

    left = adata.obs.copy()
    join_col = "__obs_join_key__"
    if use_obs_names:
        left[join_col] = adata.obs_names.astype(str)
    else:
        left[join_col] = left[obs_key].astype(str)

    meta = meta.copy()
    meta[table_key] = meta[table_key].astype(str)

    meta = meta.drop_duplicates(subset=[table_key]).copy()
    merged = left.merge(
        meta,
        left_on=join_col,
        right_on=table_key,
        how="left",
        sort=False,
    )
    merged.index = adata.obs_names
    merged["metadata_matched"] = merged[table_key].notna().to_numpy()

    drop_cols = [join_col]
    if table_key in merged.columns and (use_obs_names or table_key != obs_key):
        drop_cols.append(table_key)
    adata.obs = merged.drop(columns=[c for c in drop_cols if c in merged.columns])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    mcfg = cfg.get("metadata", {})

    if not bool(mcfg.get("enabled", False)):
        _passthrough(args.inp, args.out)
        return

    adata = sc.read_h5ad(args.inp)
    adata.obs_names_make_unique()

    if bool(mcfg.get("parse_obs_names", True)):
        _add_ids_from_obs_names(adata, mcfg)

    _merge_external_metadata(adata, mcfg)

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
