from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .utils import ensure_dir, load_config

ALIASES = {
    "cell_id": {"cellid", "cell", "cellbarcode", "barcode", "unnamed0"},
    "biological_replicate_id": {"biologicalreplicateid", "tubeid"},
    "tube_id": {"tubeid"},
    "donor_id": {"donorid", "donor", "subjectid", "participantid"},
    "age": {"age", "donorage", "ageyears"},
    "sex": {"sex", "gender"},
    "batch": {"batch", "batchid"},
    "sample_id": {"sampleid", "filename", "filepath"},
    "file_name": {"filename", "filepath"},
}


def _canonical_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _resolve_column(columns: list[str], requested: str) -> str | None:
    if requested in columns:
        return requested
    accepted = ALIASES.get(requested, {_canonical_name(requested)})
    for column in columns:
        if _canonical_name(column) in accepted:
            return column
    return None


def _clean_string(values: pd.Series) -> pd.Series:
    cleaned = values.astype("string").str.strip()
    return cleaned.mask(cleaned.eq(""))


def _replicate_config(cfg: dict[str, Any]) -> dict[str, Any]:
    section = cfg.get("biological_replicates", {})
    if not isinstance(section, dict):
        raise TypeError("biological_replicates must be a mapping.")
    return section


def canonical_replicate_column(cfg: dict[str, Any]) -> str:
    return str(
        _replicate_config(cfg).get(
            "canonical_col",
            "biological_replicate_id",
        )
    )


def attach_biological_replicates(adata, cfg: dict[str, Any]) -> dict[str, Any]:
    """Attach a configured canonical replicate ID to AnnData.obs in memory."""
    section = _replicate_config(cfg)
    if not section:
        return {"enabled": False}

    canonical_col = canonical_replicate_column(cfg)
    if canonical_col in adata.obs.columns:
        values = _clean_string(adata.obs[canonical_col])
        missing = int(values.isna().sum())
        if bool(section.get("strict", True)) and missing:
            raise ValueError(f"{canonical_col} has {missing:,} missing values in AnnData.obs.")
        adata.obs[canonical_col] = values.to_numpy()
        return {
            "enabled": True,
            "source": "adata.obs",
            "canonical_col": canonical_col,
            "matched_cells": int(values.notna().sum()),
            "missing_cells": missing,
            "n_replicates": int(values.nunique(dropna=True)),
        }

    mapping_value = section.get("mapping_path")
    if not mapping_value:
        raise FileNotFoundError(
            f"{canonical_col!r} is absent from AnnData.obs and "
            "biological_replicates.mapping_path is not configured."
        )
    mapping_path = Path(str(mapping_value))
    if not mapping_path.exists():
        raise FileNotFoundError(f"Biological-replicate mapping not found: {mapping_path}")
    print(
        f"[biological_replicates] loading cell mapping from {mapping_path}",
        flush=True,
    )

    mapping_key = str(section.get("mapping_join_key", "cell_id"))
    header = list(pd.read_csv(mapping_path, nrows=0).columns)
    key_col = _resolve_column(header, mapping_key)
    replicate_col = _resolve_column(header, canonical_col)
    if key_col is None or replicate_col is None:
        raise KeyError(
            "Could not resolve mapping columns "
            f"{mapping_key!r} and {canonical_col!r} in {mapping_path}."
        )

    mapping = pd.read_csv(
        mapping_path,
        usecols=[key_col, replicate_col],
        dtype={key_col: "string", replicate_col: "string"},
    )
    mapping[key_col] = _clean_string(mapping[key_col])
    mapping[replicate_col] = _clean_string(mapping[replicate_col])
    if bool(mapping[key_col].duplicated().any()):
        n_duplicates = int(mapping[key_col].duplicated(keep=False).sum())
        raise ValueError(f"Mapping contains {n_duplicates:,} rows with duplicate cell IDs.")

    lookup = mapping.set_index(key_col)[replicate_col]
    obs_keys = pd.Series(
        adata.obs_names.astype(str),
        index=adata.obs_names,
        dtype="string",
    )
    values = obs_keys.map(lookup)
    missing = int(values.isna().sum())
    if bool(section.get("strict", True)) and missing:
        raise ValueError(
            f"Biological-replicate mapping missed {missing:,} of {adata.n_obs:,} AnnData cells."
        )
    adata.obs[canonical_col] = values.to_numpy()
    report = {
        "enabled": True,
        "source": "mapping",
        "mapping_path": str(mapping_path),
        "mapping_join_key": key_col,
        "canonical_col": canonical_col,
        "matched_cells": int(values.notna().sum()),
        "missing_cells": missing,
        "n_replicates": int(values.nunique(dropna=True)),
    }
    print(
        "[biological_replicates] "
        f"attached {report['matched_cells']:,} cells to "
        f"{report['n_replicates']:,} replicates",
        flush=True,
    )
    return report


def add_canonical_replicate_from_obs(
    adata,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    """Create the canonical field after metadata integration on future runs."""
    section = _replicate_config(cfg)
    if not section:
        return {"enabled": False}
    canonical_col = canonical_replicate_column(cfg)
    source_col = str(section.get("source_col", "tube_id"))
    if source_col not in adata.obs.columns:
        raise KeyError(
            f"Cannot create {canonical_col}: source column "
            f"{source_col!r} is absent from AnnData.obs."
        )
    values = _clean_string(adata.obs[source_col])
    missing = int(values.isna().sum())
    if bool(section.get("strict", True)) and missing:
        raise ValueError(
            f"Cannot create {canonical_col}: {source_col} has {missing:,} missing values."
        )
    adata.obs[canonical_col] = values.to_numpy()
    return {
        "enabled": True,
        "source": f"adata.obs[{source_col}]",
        "canonical_col": canonical_col,
        "matched_cells": int(values.notna().sum()),
        "missing_cells": missing,
        "n_replicates": int(values.nunique(dropna=True)),
    }


def _conflict_count(values: pd.Series) -> int:
    return int(values.dropna().astype(str).nunique())


def _build_mapping(
    obs: pd.DataFrame,
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    section = _replicate_config(cfg)
    table_value = section.get("source_table_path")
    if not table_value:
        raise FileNotFoundError("biological_replicates.source_table_path is required.")
    table_path = Path(str(table_value))
    if not table_path.exists():
        raise FileNotFoundError(f"Source metadata table not found: {table_path}")

    sep = str(section.get("source_table_sep", ","))
    requested_key = str(section.get("source_table_join_key", "cell_id"))
    requested_source = str(section.get("source_col", "tube_id"))
    canonical_col = canonical_replicate_column(cfg)
    header = list(pd.read_csv(table_path, sep=sep, nrows=0).columns)
    key_col = _resolve_column(header, requested_key)
    source_col = _resolve_column(header, requested_source)
    if key_col is None or source_col is None:
        raise KeyError(
            f"Could not resolve source metadata columns {requested_key!r} and {requested_source!r}."
        )

    requested_metadata = {
        "donor_id": str(section.get("donor_label_col", "donor_id")),
        "age": str(section.get("age_col", "age")),
        "sex": str(section.get("sex_col", "sex")),
        "batch": str(section.get("batch_col", "batch")),
        "sample_id": str(section.get("sample_col", "file_name")),
    }
    resolved_metadata = {
        output: _resolve_column(header, requested)
        for output, requested in requested_metadata.items()
    }
    missing_metadata = [
        requested_metadata[name] for name, column in resolved_metadata.items() if column is None
    ]
    if missing_metadata:
        raise KeyError(
            "Could not resolve required source metadata columns: " + ", ".join(missing_metadata)
        )

    usecols = [key_col, source_col] + [str(column) for column in resolved_metadata.values()]
    source = pd.read_csv(
        table_path,
        sep=sep,
        usecols=list(dict.fromkeys(usecols)),
        low_memory=False,
    )
    source = source.rename(
        columns={
            key_col: "cell_id",
            source_col: canonical_col,
            **{str(column): output for output, column in resolved_metadata.items()},
        }
    )
    source["cell_id"] = _clean_string(source["cell_id"])
    source[canonical_col] = _clean_string(source[canonical_col])
    for column in ["donor_id", "sex", "batch", "sample_id"]:
        source[column] = _clean_string(source[column])
    source["age"] = pd.to_numeric(source["age"], errors="coerce")

    if source["cell_id"].isna().any():
        raise ValueError("Source metadata contains missing cell IDs.")
    if bool(source["cell_id"].duplicated().any()):
        n_duplicates = int(source["cell_id"].duplicated(keep=False).sum())
        raise ValueError(f"Source metadata contains {n_duplicates:,} duplicate cell-ID rows.")
    if not obs.index.is_unique:
        raise ValueError("AnnData obs_names are not unique.")

    source = source.set_index("cell_id")
    mapped = source.reindex(obs.index.astype(str))
    mapped.index = obs.index
    missing_cells = int(mapped[canonical_col].isna().sum())
    if bool(section.get("strict", True)) and missing_cells:
        raise ValueError(f"Source metadata missed {missing_cells:,} of {len(obs):,} cells.")

    mapping = mapped.reset_index(names="cell_id")
    replicate_audit = (
        mapping.groupby(canonical_col, observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            donor_label_nunique=("donor_id", _conflict_count),
            age_nunique=("age", _conflict_count),
            age_min=("age", "min"),
            age_max=("age", "max"),
            sex_nunique=("sex", _conflict_count),
            batch_nunique=("batch", _conflict_count),
            sample_nunique=("sample_id", _conflict_count),
        )
        .reset_index()
        .sort_values(canonical_col)
    )
    conflict_columns = [
        "donor_label_nunique",
        "age_nunique",
        "sex_nunique",
        "batch_nunique",
    ]
    conflict_counts = {
        column: int((replicate_audit[column] > 1).sum()) for column in conflict_columns
    }

    obs_consistency: dict[str, dict[str, int]] = {}
    for column in ["donor_id", "age", "sex", "batch", "sample_id"]:
        if column not in obs.columns:
            continue
        left = obs[column]
        right = mapped[column]
        if column == "age":
            left_values = pd.to_numeric(left, errors="coerce")
            right_values = pd.to_numeric(right, errors="coerce")
            comparable = left_values.notna() & right_values.notna()
            mismatches = comparable & ~np.isclose(
                left_values.fillna(0).to_numpy(dtype=float),
                right_values.fillna(0).to_numpy(dtype=float),
            )
        else:
            left_values = _clean_string(left)
            right_values = _clean_string(right)
            comparable = left_values.notna() & right_values.notna()
            mismatches = comparable & left_values.ne(right_values)
        obs_consistency[column] = {
            "comparable_cells": int(comparable.sum()),
            "mismatched_cells": int(np.asarray(mismatches).sum()),
        }

    report = {
        "source_table_path": str(table_path),
        "source_table_join_key": key_col,
        "source_replicate_column": source_col,
        "canonical_replicate_column": canonical_col,
        "n_obs": int(len(obs)),
        "mapped_cells": int(mapping[canonical_col].notna().sum()),
        "missing_cells": missing_cells,
        "unique_cell_ids": int(mapping["cell_id"].nunique()),
        "unique_biological_replicates": int(mapping[canonical_col].nunique(dropna=True)),
        "replicate_conflict_counts": conflict_counts,
        "obs_consistency": obs_consistency,
        "passed": (
            missing_cells == 0
            and int(mapping["cell_id"].nunique()) == len(obs)
            and all(value == 0 for value in conflict_counts.values())
        ),
    }
    return mapping, replicate_audit, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--mapping-out", required=True)
    parser.add_argument("--replicate-metadata-out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    adata = ad.read_h5ad(args.inp, backed="r")
    try:
        obs = adata.obs.copy()
    finally:
        adata.file.close()

    mapping, replicate_metadata, report = _build_mapping(obs, cfg)
    if bool(_replicate_config(cfg).get("strict", True)) and not report["passed"]:
        raise ValueError("Biological-replicate audit failed; inspect reported conflicts.")

    mapping_path = Path(args.mapping_out)
    replicate_path = Path(args.replicate_metadata_out)
    audit_path = Path(args.audit_out)
    ensure_dir(mapping_path.parent)
    ensure_dir(replicate_path.parent)
    ensure_dir(audit_path.parent)
    mapping.to_csv(mapping_path, index=False)
    replicate_metadata.to_csv(replicate_path, index=False)
    with audit_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        "[biological_replicates] "
        f"mapped {report['mapped_cells']:,}/{report['n_obs']:,} cells to "
        f"{report['unique_biological_replicates']:,} replicates; "
        f"passed={report['passed']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
