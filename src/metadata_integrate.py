from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import re

import numpy as np
import pandas as pd

from .biological_replicates import add_canonical_replicate_from_obs
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


def _present_mask(values: pd.Series) -> pd.Series:
    s = values.astype(str).str.strip()
    return values.notna() & s.ne("") & s.str.lower().ne("nan")


def _deduplicate_metadata(meta: pd.DataFrame, key: str) -> pd.DataFrame:
    """Collapse exact duplicates and reject ambiguous metadata keys."""
    if not _present_mask(meta[key]).all():
        n_missing = int((~_present_mask(meta[key])).sum())
        raise ValueError(
            f"metadata.table_join_key '{key}' contains {n_missing} missing values"
        )

    duplicated = meta[key].duplicated(keep=False)
    if not duplicated.any():
        return meta

    duplicate_rows = meta.loc[duplicated]
    conflicting_keys = []
    for value, group in duplicate_rows.groupby(key, sort=False, dropna=False):
        if group.drop(columns=[key]).drop_duplicates().shape[0] > 1:
            conflicting_keys.append(str(value))

    if conflicting_keys:
        preview = ", ".join(conflicting_keys[:5])
        suffix = "" if len(conflicting_keys) <= 5 else ", ..."
        raise ValueError(
            f"metadata.table_join_key '{key}' has conflicting duplicate "
            f"values for {len(conflicting_keys)} keys: {preview}{suffix}"
        )

    return meta.drop_duplicates().copy()


def _coalesce_prefer_right(merged: pd.DataFrame, base_col: str, left_tag: str, right_tag: str) -> pd.DataFrame:
    left_col = f"{base_col}_x"
    right_col = f"{base_col}_y"
    if left_col in merged.columns and right_col in merged.columns:
        left_vals = merged[left_col]
        right_vals = merged[right_col]
        right_present = _present_mask(right_vals)
        left_present = _present_mask(left_vals)
        merged[base_col] = right_vals.where(right_present, left_vals)
        merged[f"{base_col}_source"] = np.where(right_present, right_tag, np.where(left_present, left_tag, "missing"))
        both = right_present & left_present
        mismatch = both & (right_vals.astype(str) != left_vals.astype(str))
        merged[f"{base_col}_conflict"] = mismatch.to_numpy()
        n_conflicts = int(mismatch.sum())
        if n_conflicts > 0:
            print(f"[metadata] {base_col}: found {n_conflicts} row-wise conflicts between {left_col} and {right_col}.")
        return merged

    if right_col in merged.columns and base_col not in merged.columns:
        vals = merged[right_col]
        merged[base_col] = vals
        merged[f"{base_col}_source"] = np.where(_present_mask(vals), right_tag, "missing")
        merged[f"{base_col}_conflict"] = False
        return merged

    if left_col in merged.columns and base_col not in merged.columns:
        vals = merged[left_col]
        merged[base_col] = vals
        merged[f"{base_col}_source"] = np.where(_present_mask(vals), left_tag, "missing")
        merged[f"{base_col}_conflict"] = False
        return merged

    if base_col in merged.columns:
        merged[f"{base_col}_source"] = np.where(_present_mask(merged[base_col]), left_tag, "missing")
        merged[f"{base_col}_conflict"] = False
    return merged


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
    order_col = "__obs_row_order__"
    if use_obs_names:
        left[join_col] = adata.obs_names.astype(str)
    else:
        if not _present_mask(left[obs_key]).all():
            n_missing = int((~_present_mask(left[obs_key])).sum())
            raise ValueError(
                f"metadata.obs_join_key '{obs_key}' contains {n_missing} missing values"
            )
        left[join_col] = left[obs_key].astype(str)
    left[order_col] = np.arange(left.shape[0], dtype=np.int64)

    meta = meta.copy()
    meta[table_key] = meta[table_key].astype(str)

    meta = _deduplicate_metadata(meta, table_key)
    merged = left.merge(
        meta,
        left_on=join_col,
        right_on=table_key,
        how="left",
        sort=False,
        validate="many_to_one",
    )
    if merged.shape[0] != left.shape[0]:
        raise RuntimeError(
            "Metadata merge changed the number of observations: "
            f"{left.shape[0]} before, {merged.shape[0]} after"
        )
    merged = merged.sort_values(order_col, kind="stable")
    if not np.array_equal(merged[order_col].to_numpy(), np.arange(left.shape[0])):
        raise RuntimeError("Metadata merge did not preserve observation order")
    merged.index = adata.obs_names
    merged["metadata_matched"] = merged[table_key].notna().to_numpy()
    n_unmatched = int((~merged["metadata_matched"]).sum())
    if n_unmatched and bool(cfg.get("require_all_matches", True)):
        raise ValueError(
            f"Metadata join left {n_unmatched} of {left.shape[0]} observations unmatched"
        )
    if n_unmatched:
        print(
            f"[metadata] Explicitly allowing {n_unmatched} unmatched observations"
        )
    # Canonical donor_id: prefer metadata table (donor_id_y), fallback parsed (donor_id_x).
    field_donor = str(cfg.get("donor_field", "donor_id"))
    merged = _coalesce_prefer_right(
        merged,
        base_col=field_donor,
        left_tag="parsed",
        right_tag="metadata_table",
    )
    # Canonical batch: prefer metadata table (batch_y), fallback parsed (batch_x).
    field_batch = str(cfg.get("batch_field", "batch"))
    merged = _coalesce_prefer_right(
        merged,
        base_col=field_batch,
        left_tag="parsed",
        right_tag="metadata_table",
    )

    drop_cols = [join_col, order_col]
    if table_key in merged.columns and (use_obs_names or table_key != obs_key):
        drop_cols.append(table_key)
    adata.obs = merged.drop(columns=[c for c in drop_cols if c in merged.columns])


def main() -> None:
    import scanpy as sc

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
    add_canonical_replicate_from_obs(adata, cfg)

    adata.write_h5ad(args.out)


if __name__ == "__main__":
    main()
