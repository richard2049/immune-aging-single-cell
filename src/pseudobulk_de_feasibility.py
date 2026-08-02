from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from .utils import ensure_dir, load_config


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize value of type {type(value).__name__}")


def _write_json(data: dict[str, Any], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True, default=_json_default)
        handle.write("\n")


def _as_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _resolve_matrix_path(
    h5: h5py.File,
    configured_layer: str | None,
) -> tuple[str, str]:
    if configured_layer in (None, "", "X", "null", "None"):
        return "X", "adata.X"

    if configured_layer == "auto":
        if "layers" in h5 and "counts" in h5["layers"]:
            return "layers/counts", 'adata.layers["counts"]'
        return "X", "adata.X (fallback; counts layer absent)"

    path = f"layers/{configured_layer}"
    if path not in h5:
        available = sorted(h5.get("layers", {}).keys())
        raise ValueError(
            f"Configured count layer {configured_layer!r} is absent. "
            f"Available layers: {available}"
        )
    return path, f'adata.layers["{configured_layer}"]'


def _scan_csr_counts(
    h5ad_path: Path,
    configured_layer: str | None,
    rows_per_chunk: int,
    progress_every_rows: int,
    integer_tolerance: float,
) -> tuple[dict[str, Any], np.ndarray]:
    started = time.perf_counter()
    with h5py.File(h5ad_path, "r") as h5:
        matrix_path, matrix_label = _resolve_matrix_path(h5, configured_layer)
        matrix = h5[matrix_path]
        encoding = _as_text(matrix.attrs.get("encoding-type", ""))
        if encoding != "csr_matrix":
            raise ValueError(
                f"{matrix_label} uses {encoding!r}; the feasibility audit "
                "requires row-oriented CSR storage for bounded-memory scanning."
            )

        shape = tuple(int(v) for v in matrix.attrs["shape"])
        if len(shape) != 2:
            raise ValueError(f"Expected a two-dimensional matrix, found shape {shape}.")
        n_obs, n_vars = shape

        data = matrix["data"]
        indices = matrix["indices"]
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)

        structural_errors: list[str] = []
        if indptr.size != n_obs + 1:
            structural_errors.append(
                f"indptr length {indptr.size} does not equal n_obs + 1 ({n_obs + 1})"
            )
        if indptr.size and int(indptr[0]) != 0:
            structural_errors.append(f"indptr starts at {int(indptr[0])}, not 0")
        if indptr.size and int(indptr[-1]) != int(data.size):
            structural_errors.append(
                f"indptr ends at {int(indptr[-1])}, but data has {int(data.size)} values"
            )
        if int(indices.size) != int(data.size):
            structural_errors.append(
                f"indices has {int(indices.size)} entries, but data has {int(data.size)}"
            )
        if indptr.size > 1 and np.any(np.diff(indptr) < 0):
            structural_errors.append("indptr is not monotonically non-decreasing")
        if structural_errors:
            raise ValueError("Invalid CSR structure: " + "; ".join(structural_errors))

        row_library_sizes = np.zeros(n_obs, dtype=np.float64)
        nonfinite_values = 0
        negative_values = 0
        fractional_values = 0
        out_of_range_indices = 0
        total_sum = 0.0
        min_value = np.inf
        max_value = -np.inf
        max_values_in_chunk = 0
        next_progress = max(progress_every_rows, 1)

        for row_start in range(0, n_obs, rows_per_chunk):
            row_end = min(row_start + rows_per_chunk, n_obs)
            value_start = int(indptr[row_start])
            value_end = int(indptr[row_end])
            values = np.asarray(data[value_start:value_end])
            gene_indices = np.asarray(indices[value_start:value_end])
            max_values_in_chunk = max(max_values_in_chunk, int(values.size))

            if values.size:
                finite = np.isfinite(values)
                nonfinite_values += int(np.count_nonzero(~finite))
                finite_values = values[finite]
                if finite_values.size:
                    negative_values += int(np.count_nonzero(finite_values < 0))
                    fractional_values += int(
                        np.count_nonzero(
                            np.abs(finite_values - np.rint(finite_values))
                            > integer_tolerance
                        )
                    )
                    min_value = min(min_value, float(np.min(finite_values)))
                    max_value = max(max_value, float(np.max(finite_values)))
                    total_sum += float(np.sum(finite_values, dtype=np.float64))

                out_of_range_indices += int(
                    np.count_nonzero((gene_indices < 0) | (gene_indices >= n_vars))
                )

                cumulative = np.empty(values.size + 1, dtype=np.float64)
                cumulative[0] = 0.0
                np.cumsum(values, dtype=np.float64, out=cumulative[1:])
                local_ptr = indptr[row_start : row_end + 1] - value_start
                row_library_sizes[row_start:row_end] = (
                    cumulative[local_ptr[1:]] - cumulative[local_ptr[:-1]]
                )

            if row_end >= next_progress or row_end == n_obs:
                elapsed = time.perf_counter() - started
                print(
                    "[pseudobulk_feasibility] "
                    f"scanned {row_end:,}/{n_obs:,} rows "
                    f"({100.0 * row_end / n_obs:.1f}%) in {elapsed:.1f}s",
                    flush=True,
                )
                while next_progress <= row_end:
                    next_progress += max(progress_every_rows, 1)

        raw_count_checks_passed = (
            nonfinite_values == 0
            and negative_values == 0
            and fractional_values == 0
            and out_of_range_indices == 0
        )
        matrix_storage_bytes = (
            int(data.size) * int(data.dtype.itemsize)
            + int(indices.size) * int(indices.dtype.itemsize)
            + int(indptr.size) * int(matrix["indptr"].dtype.itemsize)
        )
        estimated_scan_peak_bytes = (
            max_values_in_chunk
            * (
                int(data.dtype.itemsize)
                + int(indices.dtype.itemsize)
                + np.dtype(np.float64).itemsize
            )
            + row_library_sizes.nbytes
            + indptr.nbytes
        )

        report: dict[str, Any] = {
            "matrix_path": matrix_path,
            "matrix_label": matrix_label,
            "encoding_type": encoding,
            "shape": [n_obs, n_vars],
            "data_dtype": str(data.dtype),
            "index_dtype": str(indices.dtype),
            "indptr_dtype": str(matrix["indptr"].dtype),
            "nonzero_values": int(data.size),
            "matrix_storage_bytes": matrix_storage_bytes,
            "complete_value_and_index_scan": True,
            "integer_tolerance": integer_tolerance,
            "nonfinite_values": nonfinite_values,
            "negative_values": negative_values,
            "fractional_values": fractional_values,
            "out_of_range_indices": out_of_range_indices,
            "minimum_value": None if not np.isfinite(min_value) else min_value,
            "maximum_value": None if not np.isfinite(max_value) else max_value,
            "total_count_sum": total_sum,
            "zero_library_cells": int(np.count_nonzero(row_library_sizes == 0)),
            "raw_count_checks_passed": raw_count_checks_passed,
            "rows_per_chunk": rows_per_chunk,
            "maximum_values_in_chunk": max_values_in_chunk,
            "estimated_scan_peak_bytes": int(estimated_scan_peak_bytes),
            "elapsed_seconds": time.perf_counter() - started,
        }
        return report, row_library_sizes


def _resolve_obs_column(
    columns: list[str],
    configured: str | None,
    fallbacks: list[str],
    required: bool,
) -> str | None:
    if configured and configured in columns:
        return configured
    for candidate in fallbacks:
        if candidate in columns:
            return candidate
    if required:
        requested = configured or ", ".join(fallbacks)
        raise ValueError(f"Could not resolve required obs column from: {requested}")
    return None


def _clean_string_column(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip()
    return out.mask(out.eq(""))


def _canon_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def _resolve_external_column(
    columns: list[str],
    requested: str,
) -> str | None:
    aliases = {
        "cellid": {"cellid", "cell", "cellbarcode", "barcode", "unnamed0"},
        "tubeid": {"tubeid"},
    }
    requested_canonical = _canon_column_name(requested)
    accepted = aliases.get(requested_canonical, {requested_canonical})
    for column in columns:
        if _canon_column_name(column) in accepted:
            return column
    return None


def _load_external_replicate_ids(
    obs_index: pd.Index,
    cfg: dict[str, Any],
    replicate_col: str,
) -> tuple[pd.Series, dict[str, Any]]:
    pcfg = cfg.get("pseudobulk_de_feasibility", {})
    mcfg = cfg.get("metadata", {})
    table_value = pcfg.get("replicate_table_path") or mcfg.get("table_path")
    if not table_value:
        raise FileNotFoundError(
            f"Replicate column {replicate_col!r} is absent from obs and no "
            "replicate metadata table is configured."
        )
    table_path = Path(table_value)
    if not table_path.exists():
        raise FileNotFoundError(
            f"Replicate metadata table not found: {table_path}"
        )

    sep = str(pcfg.get("replicate_table_sep", mcfg.get("table_sep", ",")))
    requested_key = str(
        pcfg.get("replicate_table_join_key")
        or mcfg.get("table_join_key", "cell_id")
    )
    header = list(pd.read_csv(table_path, sep=sep, nrows=0).columns)
    table_key = _resolve_external_column(header, requested_key)
    table_replicate_col = _resolve_external_column(header, replicate_col)
    if table_key is None or table_replicate_col is None:
        raise KeyError(
            "Could not resolve replicate metadata columns. "
            f"Requested key={requested_key!r}, replicate={replicate_col!r}; "
            f"available columns={header}"
        )

    metadata = pd.read_csv(
        table_path,
        sep=sep,
        usecols=[table_key, table_replicate_col],
        dtype={table_key: "string", table_replicate_col: "string"},
    )
    metadata[table_key] = _clean_string_column(metadata[table_key])
    metadata[table_replicate_col] = _clean_string_column(
        metadata[table_replicate_col]
    )
    metadata = metadata.dropna(subset=[table_key])

    conflicts = (
        metadata.groupby(table_key, observed=True)[table_replicate_col]
        .nunique(dropna=True)
        .gt(1)
    )
    if bool(conflicts.any()):
        raise ValueError(
            f"{int(conflicts.sum())} cell IDs map to multiple replicate IDs "
            f"in {table_path}."
        )

    lookup = (
        metadata.drop_duplicates(subset=[table_key])
        .set_index(table_key)[table_replicate_col]
    )
    replicate_ids = pd.Series(
        obs_index.astype(str),
        index=obs_index,
        dtype="string",
    ).map(lookup)
    report = {
        "source": "external_metadata",
        "table_path": str(table_path),
        "table_join_key": table_key,
        "table_replicate_col": table_replicate_col,
        "matched_cells": int(replicate_ids.notna().sum()),
        "missing_cells": int(replicate_ids.isna().sum()),
    }
    return replicate_ids, report


def _load_obs(
    h5ad_path: Path,
    config: dict[str, Any],
    expected_n_obs: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        columns = list(adata.obs.columns)
        pcfg = config.get("pseudobulk_de_feasibility", {})
        age_col = _resolve_obs_column(
            columns, pcfg.get("age_col"), ["age"], required=True
        )
        donor_col = _resolve_obs_column(
            columns,
            pcfg.get("donor_col"),
            ["donor_id", "donor", "sample_id"],
            required=True,
        )
        replicate_col = str(pcfg.get("replicate_col", donor_col))
        replicate_in_obs = replicate_col in columns and replicate_col != donor_col
        celltype_col = _resolve_obs_column(
            columns,
            pcfg.get("celltype_col"),
            ["cell_type", "majority_voting", "predicted_labels", "leiden"],
            required=True,
        )
        sample_col = _resolve_obs_column(
            columns,
            pcfg.get("sample_col"),
            ["sample_id"],
            required=False,
        )
        covariate_cols = [
            str(col)
            for col in pcfg.get("covariate_cols", ["sex", "batch"])
            if str(col) in columns
        ]

        selected = [age_col, donor_col, celltype_col]
        if replicate_in_obs and replicate_col not in selected:
            selected.append(replicate_col)
        if sample_col:
            selected.append(sample_col)
        selected.extend(col for col in covariate_cols if col not in selected)
        obs = adata.obs[selected].copy()
    finally:
        adata.file.close()

    if len(obs) != expected_n_obs:
        raise ValueError(
            f"obs has {len(obs):,} rows, but the count matrix has "
            f"{expected_n_obs:,} rows."
        )

    rename = {age_col: "age", donor_col: "donor_label", celltype_col: "cell_type"}
    if replicate_in_obs:
        rename[replicate_col] = "replicate_id"
    if sample_col:
        rename[sample_col] = "sample_id"
    obs = obs.rename(columns=rename)
    obs["age"] = pd.to_numeric(obs["age"], errors="coerce")
    obs["donor_label"] = _clean_string_column(obs["donor_label"])
    obs["cell_type"] = _clean_string_column(obs["cell_type"])
    if replicate_in_obs:
        obs["replicate_id"] = _clean_string_column(obs["replicate_id"])
        replicate_source = {
            "source": "obs",
            "obs_column": replicate_col,
            "matched_cells": int(obs["replicate_id"].notna().sum()),
            "missing_cells": int(obs["replicate_id"].isna().sum()),
        }
    elif replicate_col == donor_col:
        obs["replicate_id"] = obs["donor_label"]
        replicate_source = {
            "source": "obs_donor_label",
            "obs_column": donor_col,
            "matched_cells": int(obs["replicate_id"].notna().sum()),
            "missing_cells": int(obs["replicate_id"].isna().sum()),
        }
    else:
        replicate_ids, replicate_source = _load_external_replicate_ids(
            obs_index=obs.index,
            cfg=config,
            replicate_col=replicate_col,
        )
        obs["replicate_id"] = replicate_ids

    if "sample_id" in obs:
        obs["sample_id"] = _clean_string_column(obs["sample_id"])
        sample_source = sample_col
    else:
        obs["sample_id"] = obs["replicate_id"]
        sample_source = "replicate_id fallback"

    for col in covariate_cols:
        obs[col] = _clean_string_column(obs[col])

    resolved = {
        "age_col": age_col,
        "donor_label_col": donor_col,
        "biological_replicate_col": replicate_col,
        "biological_replicate_source": replicate_source,
        "sample_col": sample_source,
        "celltype_col": celltype_col,
        "covariate_cols": covariate_cols,
    }
    return obs, resolved


def _join_unique(values: pd.Series) -> str:
    unique = sorted(values.dropna().astype(str).unique())
    return "|".join(unique)


def _audit_metadata(
    obs: pd.DataFrame,
    covariate_cols: list[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    valid_replicate = obs.dropna(subset=["replicate_id"]).copy()
    aggregations: dict[str, tuple[str, str | Any]] = {
        "n_cells": ("replicate_id", "size"),
        "donor_label_nunique": (
            "donor_label",
            lambda x: int(x.nunique(dropna=True)),
        ),
        "donor_labels": ("donor_label", _join_unique),
        "age_nunique": ("age", lambda x: int(x.nunique(dropna=True))),
        "age_min": ("age", "min"),
        "age_max": ("age", "max"),
        "sample_nunique": ("sample_id", lambda x: int(x.nunique(dropna=True))),
        "sample_ids": ("sample_id", _join_unique),
    }
    for col in covariate_cols:
        aggregations[f"{col}_nunique"] = (
            col,
            lambda x: int(x.nunique(dropna=True)),
        )
        aggregations[f"{col}_values"] = (col, _join_unique)

    replicate_audit = (
        valid_replicate.groupby("replicate_id", observed=True)
        .agg(**aggregations)
        .reset_index()
        .sort_values("replicate_id")
    )

    pairs = obs[["replicate_id", "sample_id"]].dropna().drop_duplicates()
    samples_per_replicate = pairs.groupby("replicate_id", observed=True)[
        "sample_id"
    ].nunique()
    replicates_per_sample = pairs.groupby("sample_id", observed=True)[
        "replicate_id"
    ].nunique()
    if pairs.empty:
        relationship = "unavailable"
    elif (
        int(samples_per_replicate.max()) == 1
        and int(replicates_per_sample.max()) == 1
    ):
        relationship = "one_to_one"
    elif (
        int(samples_per_replicate.max()) > 1
        and int(replicates_per_sample.max()) == 1
    ):
        relationship = "multiple_samples_per_replicate"
    elif (
        int(samples_per_replicate.max()) == 1
        and int(replicates_per_sample.max()) > 1
    ):
        relationship = "multiple_replicates_per_sample"
    else:
        relationship = "many_to_many"

    conflict_columns = ["donor_label_nunique", "age_nunique"] + [
        f"{col}_nunique" for col in covariate_cols
    ]
    conflict_counts = {
        col.removesuffix("_nunique"): int((replicate_audit[col] > 1).sum())
        for col in conflict_columns
    }
    missing_counts = {
        col: int(obs[col].isna().sum())
        for col in [
            "replicate_id",
            "donor_label",
            "sample_id",
            "age",
            "cell_type",
            *covariate_cols,
        ]
        if col in obs
    }
    report = {
        "n_cells": int(len(obs)),
        "n_biological_replicates": int(
            obs["replicate_id"].nunique(dropna=True)
        ),
        "n_donor_labels": int(obs["donor_label"].nunique(dropna=True)),
        "n_samples": int(obs["sample_id"].nunique(dropna=True)),
        "n_cell_types": int(obs["cell_type"].nunique(dropna=True)),
        "missing_cell_values": missing_counts,
        "replicate_metadata_conflicts": conflict_counts,
        "replicate_sample_relationship": relationship,
        "max_samples_per_replicate": (
            int(samples_per_replicate.max())
            if not samples_per_replicate.empty
            else None
        ),
        "max_replicates_per_sample": (
            int(replicates_per_sample.max())
            if not replicates_per_sample.empty
            else None
        ),
    }
    return replicate_audit, report


def _build_replicate_celltype_support(
    obs: pd.DataFrame,
    row_library_sizes: np.ndarray,
    covariate_cols: list[str],
    min_cells: int,
) -> pd.DataFrame:
    working = obs.copy()
    working["library_size"] = row_library_sizes
    working = working.dropna(subset=["replicate_id", "age", "cell_type"])

    aggregations: dict[str, tuple[str, str | Any]] = {
        "n_cells": ("cell_type", "size"),
        "library_size": ("library_size", "sum"),
        "donor_label": ("donor_label", _join_unique),
        "age": ("age", "median"),
        "age_nunique": ("age", lambda x: int(x.nunique(dropna=True))),
        "n_samples": ("sample_id", lambda x: int(x.nunique(dropna=True))),
        "sample_ids": ("sample_id", _join_unique),
    }
    for col in covariate_cols:
        aggregations[col] = (col, _join_unique)
        aggregations[f"{col}_nunique"] = (
            col,
            lambda x: int(x.nunique(dropna=True)),
        )

    support = (
        working.groupby(["replicate_id", "cell_type"], observed=True)
        .agg(**aggregations)
        .reset_index()
    )
    support["passes_candidate_min_cells"] = support["n_cells"] >= min_cells
    return support.sort_values(["cell_type", "replicate_id"]).reset_index(
        drop=True
    )


def _build_design_diagnostic(
    profiles: pd.DataFrame,
    covariate_cols: list[str],
    min_residual_df: int,
) -> dict[str, Any]:
    required = ["age", *covariate_cols]
    complete = profiles.dropna(subset=required).copy()
    for col in covariate_cols:
        complete = complete[complete[col].astype("string").str.len().fillna(0) > 0]

    n_profiles = int(len(profiles))
    n_complete = int(len(complete))
    if n_complete == 0:
        return {
            "n_profiles": n_profiles,
            "n_complete_profiles": 0,
            "design_columns": 0,
            "design_rank": 0,
            "residual_df": 0,
            "full_rank": False,
            "age_span": np.nan,
            "age_r2_from_covariates": np.nan,
            "condition_number": np.nan,
            "estimable": False,
            "design_warning": "no_complete_profiles",
        }

    age = complete["age"].to_numpy(dtype=float)
    age_per_10y = (age - float(np.mean(age))) / 10.0
    covariate_parts: list[np.ndarray] = []
    covariate_labels: list[str] = []
    for col in covariate_cols:
        dummies = pd.get_dummies(
            complete[col].astype(str),
            prefix=col,
            drop_first=True,
            dtype=float,
        )
        if dummies.shape[1]:
            covariate_parts.append(dummies.to_numpy(dtype=float))
            covariate_labels.extend(str(name) for name in dummies.columns)

    intercept = np.ones((n_complete, 1), dtype=float)
    if covariate_parts:
        covariate_matrix = np.column_stack([intercept, *covariate_parts])
    else:
        covariate_matrix = intercept
    design = np.column_stack([covariate_matrix, age_per_10y])

    design_rank = int(np.linalg.matrix_rank(design))
    residual_df = int(n_complete - design_rank)
    full_rank = design_rank == design.shape[1]
    condition_number = float(np.linalg.cond(design))

    fitted_age, *_ = np.linalg.lstsq(covariate_matrix, age, rcond=None)
    residual_age = age - covariate_matrix @ fitted_age
    total_ss = float(np.sum((age - float(np.mean(age))) ** 2))
    residual_ss = float(np.sum(residual_age**2))
    age_r2 = 0.0 if total_ss == 0 else 1.0 - residual_ss / total_ss
    age_span = float(np.max(age) - np.min(age))

    warnings: list[str] = []
    if not full_rank:
        warnings.append("rank_deficient")
    if residual_df < min_residual_df:
        warnings.append("insufficient_residual_df")
    if total_ss == 0:
        warnings.append("no_age_variation")
    if age_r2 > 0.9:
        warnings.append("age_highly_explained_by_covariates")
    if not np.isfinite(condition_number) or condition_number > 1e8:
        warnings.append("ill_conditioned")

    estimable = (
        full_rank
        and residual_df >= min_residual_df
        and total_ss > 0
        and np.isfinite(condition_number)
    )
    return {
        "n_profiles": n_profiles,
        "n_complete_profiles": n_complete,
        "design_columns": int(design.shape[1]),
        "design_rank": design_rank,
        "residual_df": residual_df,
        "full_rank": full_rank,
        "age_span": age_span,
        "age_r2_from_covariates": float(age_r2),
        "condition_number": condition_number,
        "covariate_dummy_columns": "|".join(covariate_labels),
        "estimable": estimable,
        "design_warning": "|".join(warnings),
    }


def _summarize_celltype_support(
    support: pd.DataFrame,
    covariate_cols: list[str],
    min_cells: int,
    min_donors: int,
    min_age_span: float,
    min_residual_df: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, Any]] = []
    design_rows: list[dict[str, Any]] = []

    for cell_type, all_profiles in support.groupby("cell_type", observed=True):
        passing = all_profiles[all_profiles["n_cells"] >= min_cells].copy()
        design = _build_design_diagnostic(
            passing,
            covariate_cols=covariate_cols,
            min_residual_df=min_residual_df,
        )
        design_row = {"cell_type": str(cell_type), **design}
        design_rows.append(design_row)

        age_min = float(passing["age"].min()) if not passing.empty else np.nan
        age_max = float(passing["age"].max()) if not passing.empty else np.nan
        age_span = age_max - age_min if not passing.empty else np.nan
        exclusion_reasons: list[str] = []
        if len(passing) < min_donors:
            exclusion_reasons.append(
                "insufficient_biological_replicates_with_min_cells"
            )
        if not np.isfinite(age_span) or age_span < min_age_span:
            exclusion_reasons.append("insufficient_age_span")
        if not bool(design["estimable"]):
            exclusion_reasons.append("design_not_estimable")

        summary_rows.append(
            {
                "cell_type": str(cell_type),
                "total_cells": int(all_profiles["n_cells"].sum()),
                "replicates_total": int(
                    all_profiles["replicate_id"].nunique()
                ),
                "replicates_meeting_min_cells": int(
                    passing["replicate_id"].nunique()
                ),
                "candidate_min_cells": min_cells,
                "median_cells_per_passing_donor": (
                    float(passing["n_cells"].median())
                    if not passing.empty
                    else np.nan
                ),
                "median_library_size_per_passing_donor": (
                    float(passing["library_size"].median())
                    if not passing.empty
                    else np.nan
                ),
                "age_min": age_min,
                "age_max": age_max,
                "age_span": age_span,
                "candidate_min_donors": min_donors,
                "candidate_min_age_span": min_age_span,
                "candidate_eligible": not exclusion_reasons,
                "exclusion_reason": "|".join(exclusion_reasons),
            }
        )

    summary = pd.DataFrame(summary_rows)
    design_df = pd.DataFrame(design_rows)
    summary = summary.merge(
        design_df[
            [
                "cell_type",
                "n_complete_profiles",
                "design_columns",
                "design_rank",
                "residual_df",
                "full_rank",
                "age_r2_from_covariates",
                "condition_number",
                "estimable",
                "design_warning",
            ]
        ],
        on="cell_type",
        how="left",
    )
    summary = summary.sort_values(
        [
            "candidate_eligible",
            "replicates_meeting_min_cells",
            "total_cells",
        ],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    design_df = design_df.sort_values("cell_type").reset_index(drop=True)
    return summary, design_df


def _audit_edger_environment() -> dict[str, Any]:
    rscript = shutil.which("Rscript")
    report: dict[str, Any] = {
        "rscript_path": rscript,
        "edger_available": False,
        "edger_version": None,
    }
    if not rscript:
        return report

    try:
        result = subprocess.run(
            [
                rscript,
                "--vanilla",
                "-e",
                "cat(as.character(packageVersion('edgeR')))",
            ],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        report["error"] = str(exc)
        return report

    report["return_code"] = int(result.returncode)
    if result.returncode == 0:
        report["edger_available"] = True
        report["edger_version"] = result.stdout.strip()
    else:
        report["error"] = result.stderr.strip()
    return report


def run_feasibility_audit(
    config_path: Path,
    h5ad_path: Path,
    outdir: Path,
) -> dict[str, Any]:
    started = time.perf_counter()
    cfg = load_config(config_path)
    pcfg = cfg.get("pseudobulk_de_feasibility", {})
    ensure_dir(outdir)

    min_cells = int(pcfg.get("candidate_min_cells_per_pseudobulk", 50))
    min_donors = int(pcfg.get("candidate_min_donors_per_celltype", 12))
    min_age_span = float(pcfg.get("candidate_min_age_span", 12.0))
    min_residual_df = int(pcfg.get("candidate_min_residual_df", 5))
    rows_per_chunk = int(pcfg.get("scan_rows_per_chunk", 2000))
    progress_every_rows = int(pcfg.get("progress_every_rows", 100000))
    integer_tolerance = float(pcfg.get("integer_tolerance", 1e-6))
    configured_layer = pcfg.get("count_layer")

    print(
        "[pseudobulk_feasibility] "
        f"starting full count audit: {h5ad_path}",
        flush=True,
    )
    matrix_audit, row_library_sizes = _scan_csr_counts(
        h5ad_path=h5ad_path,
        configured_layer=configured_layer,
        rows_per_chunk=rows_per_chunk,
        progress_every_rows=progress_every_rows,
        integer_tolerance=integer_tolerance,
    )
    _write_json(matrix_audit, outdir / "matrix_audit.json")
    np.save(outdir / "cell_library_sizes.npy", row_library_sizes)

    print("[pseudobulk_feasibility] loading obs metadata only", flush=True)
    obs, resolved_columns = _load_obs(
        h5ad_path=h5ad_path,
        config=cfg,
        expected_n_obs=int(matrix_audit["shape"][0]),
    )
    replicate_audit, metadata_audit = _audit_metadata(
        obs, covariate_cols=resolved_columns["covariate_cols"]
    )
    replicate_audit.to_csv(
        outdir / "replicate_metadata_audit.csv", index=False
    )

    support = _build_replicate_celltype_support(
        obs=obs,
        row_library_sizes=row_library_sizes,
        covariate_cols=resolved_columns["covariate_cols"],
        min_cells=min_cells,
    )
    support.to_csv(
        outdir / "replicate_celltype_support.csv", index=False
    )
    celltype_summary, design_diagnostics = _summarize_celltype_support(
        support=support,
        covariate_cols=resolved_columns["covariate_cols"],
        min_cells=min_cells,
        min_donors=min_donors,
        min_age_span=min_age_span,
        min_residual_df=min_residual_df,
    )
    celltype_summary.to_csv(
        outdir / "celltype_support_summary.csv", index=False
    )
    design_diagnostics.to_csv(
        outdir / "design_diagnostics.csv", index=False
    )

    edger_environment = _audit_edger_environment()
    metadata_conflicts = metadata_audit["replicate_metadata_conflicts"]
    blockers: list[str] = []
    if not matrix_audit["raw_count_checks_passed"]:
        blockers.append("raw_count_matrix_validation_failed")
    if any(int(value) > 0 for value in metadata_conflicts.values()):
        blockers.append("replicate_metadata_conflicts_present")
    if metadata_audit["missing_cell_values"].get("replicate_id", 0) > 0:
        blockers.append("missing_biological_replicate_ids")
    sample_role = str(pcfg.get("sample_role", "unknown"))
    if (
        metadata_audit["replicate_sample_relationship"] != "one_to_one"
        and sample_role != "multiplexed_technical_library"
    ):
        blockers.append("replicate_sample_mapping_requires_design_decision")
    if not bool(celltype_summary["candidate_eligible"].any()):
        blockers.append("no_cell_types_pass_candidate_support_and_design")
    if not edger_environment["edger_available"]:
        blockers.append("edger_runtime_not_available")
    blockers.append("analysis_contract_requires_human_approval")

    eligible_profiles = support[
        support["cell_type"].isin(
            celltype_summary.loc[
                celltype_summary["candidate_eligible"], "cell_type"
            ]
        )
        & support["passes_candidate_min_cells"]
    ]
    n_vars = int(matrix_audit["shape"][1])
    dense_pseudobulk_bytes = int(len(eligible_profiles) * n_vars * 8)

    summary: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path),
        "input_h5ad": str(h5ad_path),
        "input_size_bytes": int(h5ad_path.stat().st_size),
        "output_directory": str(outdir),
        "resolved_columns": resolved_columns,
        "candidate_thresholds": {
            "min_cells_per_pseudobulk": min_cells,
            "min_donors_per_celltype": min_donors,
            "min_age_span": min_age_span,
            "min_residual_df": min_residual_df,
        },
        "matrix_audit": matrix_audit,
        "metadata_audit": metadata_audit,
        "celltype_support": {
            "cell_types_total": int(len(celltype_summary)),
            "cell_types_candidate_eligible": int(
                celltype_summary["candidate_eligible"].sum()
            ),
            "eligible_cell_types": celltype_summary.loc[
                celltype_summary["candidate_eligible"], "cell_type"
            ].tolist(),
            "eligible_replicate_celltype_profiles": int(
                len(eligible_profiles)
            ),
        },
        "sample_role": sample_role,
        "memory_estimate": {
            "dense_int64_pseudobulk_matrix_bytes": dense_pseudobulk_bytes,
            "planned_aggregation": (
                "sparse chunked sum; dense cell-by-gene conversion prohibited"
            ),
        },
        "edger_environment": edger_environment,
        "implementation_blockers": blockers,
        "ready_for_de_implementation": len(blockers) == 0,
        "elapsed_seconds": time.perf_counter() - started,
    }
    _write_json(summary, outdir / "feasibility_summary.json")
    print(
        "[pseudobulk_feasibility] "
        f"completed in {summary['elapsed_seconds']:.1f}s; "
        f"candidate eligible cell types: "
        f"{summary['celltype_support']['cell_types_candidate_eligible']}; "
        f"outputs: {outdir}",
        flush=True,
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Audit raw counts, donor metadata, cell-type support, and candidate "
            "design matrices before pseudobulk differential expression."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    run_feasibility_audit(
        config_path=Path(args.config),
        h5ad_path=Path(args.inp),
        outdir=Path(args.outdir),
    )


if __name__ == "__main__":
    main()
