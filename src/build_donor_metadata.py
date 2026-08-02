from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .utils import ensure_dir


def _read_table(path: Path, sep: str | None, sheet: str | None) -> pd.DataFrame:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, sheet_name=sheet or 0)
    if sep is not None:
        return pd.read_csv(path, sep=sep)
    if suffix.endswith((".tsv", ".tsv.gz", ".txt")):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _require_columns(table: pd.DataFrame, columns: list[str], source: Path) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise KeyError(f"Missing columns in {source}: {missing}")


def _deduplicate_cell_metadata(table: pd.DataFrame, cell_col: str) -> pd.DataFrame:
    if table[cell_col].isna().any():
        raise ValueError(f"Cell join column '{cell_col}' contains missing values")
    table = table.copy()
    table[cell_col] = table[cell_col].astype(str)
    duplicated = table[cell_col].duplicated(keep=False)
    conflicts = []
    for cell_id, group in table.loc[duplicated].groupby(cell_col, sort=False):
        if group.drop(columns=[cell_col]).drop_duplicates().shape[0] > 1:
            conflicts.append(str(cell_id))
    if conflicts:
        preview = ", ".join(conflicts[:5])
        raise ValueError(
            f"Conflicting metadata rows for {len(conflicts)} cell IDs: {preview}"
        )
    return table.drop_duplicates().copy()


def _single_value(values: pd.Series, field: str, replicate_id: str):
    observed = values.dropna().drop_duplicates()
    if observed.shape[0] > 1:
        preview = ", ".join(observed.astype(str).head(5))
        raise ValueError(
            f"Biological replicate '{replicate_id}' has conflicting {field} "
            f"values: {preview}"
        )
    return np.nan if observed.empty else observed.iloc[0]


def _aggregate_replicates(cells: pd.DataFrame) -> pd.DataFrame:
    key = "biological_replicate_id"
    if cells[key].isna().any() or cells[key].astype(str).str.strip().eq("").any():
        raise ValueError("Every included cell must have a biological_replicate_id")
    cells = cells.copy()
    cells[key] = cells[key].astype(str)

    rows = []
    metadata_fields = [
        field
        for field in ("donor_id", "age", "sex", "cohort")
        if field in cells.columns
    ]
    for replicate_id, group in cells.groupby(key, sort=True):
        row = {
            key: replicate_id,
            "n_cells": int(group.shape[0]),
            "n_samples": int(group["sample_id"].nunique()),
            "batches": ";".join(
                sorted({str(value) for value in group["batch"].dropna()})
            ),
        }
        for field in metadata_fields:
            row[field] = _single_value(group[field], field, replicate_id)
        rows.append(row)

    result = pd.DataFrame(rows)
    result["has_age"] = result["age"].notna() if "age" in result else False
    result["has_sex"] = result["sex"].notna() if "sex" in result else False
    return result


def _is_raw_output(path: Path) -> bool:
    parts = [part.lower() for part in path.resolve().parts]
    return any(parts[index : index + 2] == ["data", "raw"] for index in range(len(parts) - 1))


def _cells_from_h5ad(path: Path, replicate_obs_col: str | None) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    try:
        obs_names = pd.Index(adata.obs_names.astype(str))
        cells = pd.DataFrame({"cell_id": obs_names})
        cells["sample_id"] = cells["cell_id"].str.split("_", n=1).str[0]
        cells["batch"] = cells["sample_id"]
        if replicate_obs_col:
            if replicate_obs_col not in adata.obs.columns:
                raise KeyError(
                    f"Replicate column '{replicate_obs_col}' is missing from adata.obs"
                )
            cells["biological_replicate_id"] = adata.obs[
                replicate_obs_col
            ].astype(str).to_numpy()
        return cells
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def _external_cell_metadata(args: argparse.Namespace) -> pd.DataFrame:
    if not args.supp or not args.cell_col or not args.replicate_col:
        raise ValueError(
            "Without --replicate-obs-col, provide --supp, --cell-col, and "
            "--replicate-col so the biological replicate comes from explicit metadata"
        )

    selected = {
        args.cell_col: "cell_id",
        args.replicate_col: "biological_replicate_id",
    }
    optional = {
        args.sample_col: "sample_id",
        args.batch_col: "batch",
        args.donor_col: "donor_id",
        args.age_col: "age",
        args.sex_col: "sex",
        args.cohort_col: "cohort",
    }
    selected.update({source: target for source, target in optional.items() if source})

    tables = []
    for value in args.supp:
        path = Path(value)
        if not path.exists():
            raise FileNotFoundError(f"Supplementary metadata file not found: {path}")
        table = _read_table(path, sep=args.supp_sep, sheet=args.supp_sheet)
        _require_columns(table, list(selected), path)
        tables.append(table.loc[:, list(selected)].rename(columns=selected))

    combined = pd.concat(tables, ignore_index=True)
    return _deduplicate_cell_metadata(combined, "cell_id")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build replicate-level metadata from an explicit biological-replicate key."
    )
    parser.add_argument("--h5ad", required=True, help="Input AnnData used by the pipeline.")
    parser.add_argument("--out", required=True, help="Derived replicate metadata CSV.")
    parser.add_argument("--supp", action="append", default=[], help="Cell-level metadata table.")
    parser.add_argument("--supp-sep", default=None)
    parser.add_argument("--supp-sheet", default=None)
    parser.add_argument("--replicate-obs-col", default=None)
    parser.add_argument("--cell-col", default=None)
    parser.add_argument("--replicate-col", default=None)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--batch-col", default=None)
    parser.add_argument("--donor-col", default=None)
    parser.add_argument("--age-col", default=None)
    parser.add_argument("--sex-col", default=None)
    parser.add_argument("--cohort-col", default=None)
    parser.add_argument("--allow-unmatched-cells", action="store_true")
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    out_path = Path(args.out)
    if _is_raw_output(out_path):
        raise ValueError(
            "Derived metadata must not be written under data/raw; use data/derived or results"
        )

    cells = _cells_from_h5ad(h5ad_path, args.replicate_obs_col)
    if not args.replicate_obs_col:
        metadata = _external_cell_metadata(args)
        base = cells.drop(columns=["sample_id", "batch"])
        cells = base.merge(
            metadata,
            on="cell_id",
            how="left",
            sort=False,
            validate="one_to_one",
        )
        if "sample_id" not in cells:
            cells["sample_id"] = cells["cell_id"].str.split("_", n=1).str[0]
        if "batch" not in cells:
            cells["batch"] = cells["sample_id"]
        unmatched = int(cells["biological_replicate_id"].isna().sum())
        if unmatched and not args.allow_unmatched_cells:
            raise ValueError(
                f"{unmatched} AnnData cells lack external biological-replicate metadata"
            )
        cells = cells.dropna(subset=["biological_replicate_id"]).copy()

    replicate_metadata = _aggregate_replicates(cells)
    ensure_dir(out_path.parent)
    replicate_metadata.to_csv(out_path, index=False)

    print(f"Wrote replicate metadata: {out_path}")
    print(f"Biological replicates: {replicate_metadata.shape[0]}")
    if "age" in replicate_metadata:
        print(
            "Biological replicates with age: "
            f"{int(replicate_metadata['age'].notna().sum())}"
        )


if __name__ == "__main__":
    main()
