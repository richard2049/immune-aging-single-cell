from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

from .biological_replicates import attach_biological_replicates
from .utils import ensure_dir, load_config


def _clean_string(values: pd.Series) -> pd.Series:
    cleaned = values.astype("string").str.strip()
    return cleaned.mask(cleaned.eq(""))


def _require_contract(cfg: dict[str, Any]) -> dict[str, Any]:
    contract = cfg.get("pseudobulk_de", {})
    if not isinstance(contract, dict) or not contract:
        raise KeyError("pseudobulk_de configuration is required.")

    required = [
        "replicate_col",
        "celltype_col",
        "age_col",
        "sex_col",
        "batch_col",
        "primary_cell_types",
        "exploratory_cell_types",
    ]
    missing = [key for key in required if key not in contract]
    if missing:
        raise KeyError(
            "Missing pseudobulk_de configuration keys: " + ", ".join(missing)
        )

    primary = [str(value) for value in contract["primary_cell_types"]]
    exploratory = [
        str(value) for value in contract["exploratory_cell_types"]
    ]
    overlap = sorted(set(primary).intersection(exploratory))
    if overlap:
        raise ValueError(
            "Cell types cannot be both primary and exploratory: "
            + ", ".join(overlap)
        )
    if len(primary) != len(set(primary)) or len(exploratory) != len(
        set(exploratory)
    ):
        raise ValueError("Configured cell-type lists contain duplicates.")
    return contract


def _replicate_metadata(
    obs: pd.DataFrame,
    replicate_col: str,
    age_col: str,
    sex_col: str,
    batch_col: str,
) -> pd.DataFrame:
    grouped = obs.groupby(replicate_col, observed=True, sort=True)
    conflicts = grouped[[age_col, sex_col, batch_col]].nunique(dropna=False)
    bad = conflicts.gt(1).any(axis=1)
    if bool(bad.any()):
        examples = ", ".join(conflicts.index[bad].astype(str)[:5])
        raise ValueError(
            "Age, sex, or batch is inconsistent within biological replicate; "
            f"examples: {examples}"
        )

    metadata = grouped[[age_col, sex_col, batch_col]].first().reset_index()
    metadata = metadata.rename(
        columns={
            replicate_col: "biological_replicate_id",
            age_col: "age",
            sex_col: "sex",
            batch_col: "batch",
        }
    )
    metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce")
    if metadata[["age", "sex", "batch"]].isna().any().any():
        raise ValueError("Replicate-level age, sex, and batch must be complete.")
    metadata["age_decade"] = metadata["age"] / 10.0
    return metadata


def _design_preflight(metadata: pd.DataFrame) -> dict[str, Any]:
    n_profiles = int(len(metadata))
    factor_levels = {
        "sex": int(metadata["sex"].nunique(dropna=True)),
        "batch": int(metadata["batch"].nunique(dropna=True)),
    }
    age_values = metadata["age_decade"].to_numpy(dtype=float)
    age_varies = bool(n_profiles > 1 and np.ptp(age_values) > 0)

    encoded = pd.get_dummies(
        metadata[["sex", "batch"]].astype("string"),
        columns=["sex", "batch"],
        drop_first=True,
        dtype=float,
    )
    design = np.column_stack(
        [
            np.ones(n_profiles, dtype=float),
            encoded.to_numpy(dtype=float),
            age_values,
        ]
    )
    rank = int(np.linalg.matrix_rank(design)) if n_profiles else 0
    n_columns = int(design.shape[1])
    estimable = (
        n_profiles > 0
        and age_varies
        and all(value >= 2 for value in factor_levels.values())
        and rank == n_columns
    )
    return {
        "design_columns": n_columns,
        "design_rank": rank,
        "residual_df": n_profiles - rank,
        "sex_levels": factor_levels["sex"],
        "batch_levels": factor_levels["batch"],
        "full_design_estimable": estimable,
    }


def _build_profiles_and_eligibility(
    obs: pd.DataFrame,
    replicate_metadata: pd.DataFrame,
    contract: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    replicate_col = str(contract["replicate_col"])
    celltype_col = str(contract["celltype_col"])
    sample_col = str(contract.get("sample_col", ""))
    min_cells = int(contract.get("min_cells_per_pseudobulk", 50))
    min_replicates = int(contract.get("min_replicates_per_celltype", 12))
    min_age_span = float(contract.get("min_age_span_years", 12))
    min_residual_df = int(contract.get("min_residual_df", 5))

    counts = (
        obs.groupby([celltype_col, replicate_col], observed=True, sort=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(
            columns={
                celltype_col: "cell_type",
                replicate_col: "biological_replicate_id",
            }
        )
    )
    counts = counts.merge(
        replicate_metadata,
        on="biological_replicate_id",
        how="left",
        validate="many_to_one",
    )

    if sample_col and sample_col in obs.columns:
        technical_libraries = (
            obs.groupby(
                [celltype_col, replicate_col],
                observed=True,
                sort=True,
            )[sample_col]
            .nunique(dropna=True)
            .rename("n_technical_libraries")
            .reset_index()
            .rename(
                columns={
                    celltype_col: "cell_type",
                    replicate_col: "biological_replicate_id",
                }
            )
        )
        counts = counts.merge(
            technical_libraries,
            on=["cell_type", "biological_replicate_id"],
            how="left",
            validate="one_to_one",
        )
    else:
        counts["n_technical_libraries"] = pd.NA

    primary = {str(value) for value in contract["primary_cell_types"]}
    exploratory = {
        str(value) for value in contract["exploratory_cell_types"]
    }
    tier = {
        **{cell_type: "primary" for cell_type in primary},
        **{cell_type: "exploratory" for cell_type in exploratory},
    }

    eligibility_rows: list[dict[str, Any]] = []
    included_cell_types: set[str] = set()
    for cell_type, group in counts.groupby("cell_type", observed=True):
        qualifying = group.loc[group["n_cells"] >= min_cells].copy()
        age_span = (
            float(qualifying["age"].max() - qualifying["age"].min())
            if not qualifying.empty
            else np.nan
        )
        design = _design_preflight(qualifying)
        reasons: list[str] = []
        analysis_tier = tier.get(str(cell_type), "not_approved")
        if analysis_tier == "not_approved":
            reasons.append("not_approved_by_analysis_contract")
        if len(qualifying) < min_replicates:
            reasons.append("insufficient_qualifying_replicates")
        if not np.isfinite(age_span) or age_span < min_age_span:
            reasons.append("insufficient_age_span")
        if design["residual_df"] < min_residual_df:
            reasons.append("insufficient_full_model_residual_df")

        support_passes = not any(
            reason
            in {
                "insufficient_qualifying_replicates",
                "insufficient_age_span",
                "insufficient_full_model_residual_df",
            }
            for reason in reasons
        )
        aggregate = analysis_tier != "not_approved" and support_passes
        if aggregate:
            included_cell_types.add(str(cell_type))

        eligibility_rows.append(
            {
                "cell_type": str(cell_type),
                "analysis_tier": analysis_tier,
                "n_replicates_observed": int(len(group)),
                "n_replicates_qualifying": int(len(qualifying)),
                "age_span_years": age_span,
                **design,
                "aggregation_status": (
                    "included" if aggregate else "excluded"
                ),
                "exclusion_reasons": ";".join(reasons),
            }
        )

    configured_missing = sorted((primary | exploratory) - set(counts["cell_type"]))
    for cell_type in configured_missing:
        eligibility_rows.append(
            {
                "cell_type": cell_type,
                "analysis_tier": tier[cell_type],
                "n_replicates_observed": 0,
                "n_replicates_qualifying": 0,
                "age_span_years": np.nan,
                "design_columns": 0,
                "design_rank": 0,
                "residual_df": 0,
                "sex_levels": 0,
                "batch_levels": 0,
                "full_design_estimable": False,
                "aggregation_status": "excluded",
                "exclusion_reasons": "configured_cell_type_not_observed",
            }
        )

    profiles = counts.loc[
        counts["cell_type"].isin(included_cell_types)
        & counts["n_cells"].ge(min_cells)
    ].copy()
    profiles["analysis_tier"] = profiles["cell_type"].map(tier)
    profiles = profiles.sort_values(
        ["cell_type", "biological_replicate_id"],
        kind="stable",
    ).reset_index(drop=True)
    profiles.insert(
        0,
        "profile_id",
        [f"PB{index:05d}" for index in range(1, len(profiles) + 1)],
    )
    eligibility = pd.DataFrame(eligibility_rows).sort_values(
        ["analysis_tier", "cell_type"],
        kind="stable",
    )
    return profiles, eligibility


def _profile_codes(
    obs: pd.DataFrame,
    profiles: pd.DataFrame,
    replicate_col: str,
    celltype_col: str,
) -> np.ndarray:
    profile_index = pd.MultiIndex.from_frame(
        profiles[["biological_replicate_id", "cell_type"]]
    )
    cell_index = pd.MultiIndex.from_arrays(
        [
            obs[replicate_col].astype(str),
            obs[celltype_col].astype(str),
        ],
        names=["biological_replicate_id", "cell_type"],
    )
    return profile_index.get_indexer(cell_index).astype(np.int64, copy=False)


def _integer_sparse_chunk(matrix: Any) -> sp.csr_matrix:
    if not sp.issparse(matrix):
        raise TypeError(
            "Pseudobulk aggregation requires a sparse count matrix; a dense "
            "cell-by-gene matrix is not permitted."
        )
    chunk = matrix.tocsr()
    values = np.asarray(chunk.data)
    if values.size:
        if not bool(np.isfinite(values).all()):
            raise ValueError("Count matrix contains non-finite values.")
        if bool((values < 0).any()):
            raise ValueError("Count matrix contains negative values.")
        rounded = np.rint(values)
        if not bool(np.allclose(values, rounded, rtol=0, atol=1e-6)):
            raise ValueError("Count matrix contains fractional values.")
        chunk.data = rounded.astype(np.int64, copy=False)
    return chunk.astype(np.int64, copy=False)


def aggregate_sparse_counts(
    matrix: Any,
    profile_codes: np.ndarray,
    n_profiles: int,
    chunk_size: int,
) -> tuple[sp.csr_matrix, np.ndarray, int]:
    n_obs, n_genes = matrix.shape
    if len(profile_codes) != n_obs:
        raise ValueError("Profile-code vector does not match count-matrix rows.")
    if n_profiles <= 0:
        raise ValueError("No eligible pseudobulk profiles were defined.")

    aggregate = sp.csr_matrix((n_profiles, n_genes), dtype=np.int64)
    expected_library_sizes = np.zeros(n_profiles, dtype=np.int64)
    included_cells = 0

    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        codes = profile_codes[start:end]
        selected = codes >= 0
        if not bool(selected.any()):
            continue

        chunk = _integer_sparse_chunk(matrix[start:end])
        selected_chunk = chunk[selected]
        selected_codes = codes[selected]
        indicator = sp.csr_matrix(
            (
                np.ones(len(selected_codes), dtype=np.int64),
                (selected_codes, np.arange(len(selected_codes))),
            ),
            shape=(n_profiles, len(selected_codes)),
        )
        aggregate = aggregate + indicator @ selected_chunk
        row_sums = np.asarray(selected_chunk.sum(axis=1)).ravel()
        expected_library_sizes += np.bincount(
            selected_codes,
            weights=row_sums,
            minlength=n_profiles,
        ).astype(np.int64)
        included_cells += int(selected.sum())

        print(
            "[pseudobulk_aggregate] "
            f"processed {end:,}/{n_obs:,} cells; "
            f"included={included_cells:,}",
            flush=True,
        )

    aggregate.sum_duplicates()
    aggregate.eliminate_zeros()
    aggregate.sort_indices()
    observed_library_sizes = np.asarray(aggregate.sum(axis=1)).ravel().astype(
        np.int64
    )
    if not np.array_equal(observed_library_sizes, expected_library_sizes):
        raise RuntimeError(
            "Aggregated library sizes do not match selected raw-count totals."
        )
    return aggregate, observed_library_sizes, included_cells


def _write_csv_atomic(table: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    table.to_csv(temporary, index=False)
    os.replace(temporary, path)


def _write_matrix_market_gzip_atomic(
    matrix: sp.csr_matrix,
    path: Path,
) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw,
            mtime=0,
        ) as compressed:
            mmwrite(
                compressed,
                matrix,
                field="integer",
                symmetry="general",
            )
    os.replace(temporary, path)


def run(
    config_path: Path,
    input_h5ad: Path,
    matrix_out: Path,
    profiles_out: Path,
    genes_out: Path,
    eligibility_out: Path,
    audit_out: Path,
) -> dict[str, Any]:
    cfg = load_config(config_path)
    contract = _require_contract(cfg)
    chunk_size = int(contract.get("aggregation_chunk_size", 5000))
    if chunk_size <= 0:
        raise ValueError("aggregation_chunk_size must be positive.")

    adata = ad.read_h5ad(input_h5ad, backed="r")
    try:
        mapping_report = attach_biological_replicates(adata, cfg)
        obs = adata.obs.copy()
        replicate_col = str(contract["replicate_col"])
        celltype_col = str(contract["celltype_col"])
        age_col = str(contract["age_col"])
        sex_col = str(contract["sex_col"])
        batch_col = str(contract["batch_col"])
        required_columns = [
            replicate_col,
            celltype_col,
            age_col,
            sex_col,
            batch_col,
        ]
        missing = [column for column in required_columns if column not in obs]
        if missing:
            raise KeyError(
                "Annotated checkpoint lacks required columns: "
                + ", ".join(missing)
            )

        obs[replicate_col] = _clean_string(obs[replicate_col])
        obs[celltype_col] = _clean_string(obs[celltype_col])
        obs[sex_col] = _clean_string(obs[sex_col])
        obs[batch_col] = _clean_string(obs[batch_col])
        obs[age_col] = pd.to_numeric(obs[age_col], errors="coerce")
        if obs[required_columns].isna().any().any():
            raise ValueError(
                "Pseudobulk grouping and model metadata must be complete."
            )

        replicate_metadata = _replicate_metadata(
            obs,
            replicate_col=replicate_col,
            age_col=age_col,
            sex_col=sex_col,
            batch_col=batch_col,
        )
        profiles, eligibility = _build_profiles_and_eligibility(
            obs,
            replicate_metadata,
            contract,
        )
        codes = _profile_codes(
            obs,
            profiles,
            replicate_col=replicate_col,
            celltype_col=celltype_col,
        )

        count_layer = contract.get("count_layer")
        if count_layer in {None, "", "X"}:
            matrix = adata.X
            matrix_source = "adata.X"
        else:
            if str(count_layer) not in adata.layers:
                raise KeyError(f"Count layer not found: {count_layer}")
            matrix = adata.layers[str(count_layer)]
            matrix_source = f"adata.layers[{count_layer!r}]"

        aggregate, library_sizes, included_cells = aggregate_sparse_counts(
            matrix,
            profile_codes=codes,
            n_profiles=len(profiles),
            chunk_size=chunk_size,
        )
        profiles["library_size"] = library_sizes

        gene_ids = pd.Index(adata.var_names.astype(str), name="gene_id")
        if not gene_ids.is_unique:
            raise ValueError(
                "Gene identifiers must be unique before pseudobulk DE."
            )
        genes = pd.DataFrame({"gene_id": gene_ids})
    finally:
        adata.file.close()

    _write_matrix_market_gzip_atomic(aggregate, matrix_out)
    _write_csv_atomic(profiles, profiles_out)
    _write_csv_atomic(genes, genes_out)
    _write_csv_atomic(eligibility, eligibility_out)

    audit = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path),
        "input_h5ad": str(input_h5ad),
        "matrix_source": matrix_source,
        "matrix_orientation": "profiles_by_genes",
        "matrix_shape": [int(value) for value in aggregate.shape],
        "matrix_nonzero_values": int(aggregate.nnz),
        "aggregated_total_counts": int(aggregate.sum()),
        "included_cells": included_cells,
        "excluded_cells": int(len(obs) - included_cells),
        "n_profiles": int(len(profiles)),
        "n_primary_profiles": int(
            profiles["analysis_tier"].eq("primary").sum()
        ),
        "n_exploratory_profiles": int(
            profiles["analysis_tier"].eq("exploratory").sum()
        ),
        "n_genes": int(len(genes)),
        "mapping_report": mapping_report,
        "approved_contract": {
            key: contract.get(key)
            for key in [
                "min_cells_per_pseudobulk",
                "min_replicates_per_celltype",
                "min_age_span_years",
                "min_residual_df",
                "primary_formula",
                "age_effect_scale",
                "primary_cell_types",
                "exploratory_cell_types",
            ]
        },
        "count_conservation_passed": True,
        "dense_cell_by_gene_conversion": False,
        "output_files": {
            "matrix": str(matrix_out),
            "profiles": str(profiles_out),
            "genes": str(genes_out),
            "eligibility": str(eligibility_out),
        },
    }
    ensure_dir(audit_out.parent)
    temporary_audit = audit_out.with_name(audit_out.name + ".tmp")
    with temporary_audit.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary_audit, audit_out)
    print(
        "[pseudobulk_aggregate] "
        f"wrote {len(profiles):,} profiles x {len(genes):,} genes; "
        f"total_counts={int(aggregate.sum()):,}",
        flush=True,
    )
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate raw sparse counts by biological replicate and cell type "
            "under the approved pseudobulk DE contract."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--matrix-out", required=True)
    parser.add_argument("--profiles-out", required=True)
    parser.add_argument("--genes-out", required=True)
    parser.add_argument("--eligibility-out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()

    run(
        config_path=Path(args.config),
        input_h5ad=Path(args.inp),
        matrix_out=Path(args.matrix_out),
        profiles_out=Path(args.profiles_out),
        genes_out=Path(args.genes_out),
        eligibility_out=Path(args.eligibility_out),
        audit_out=Path(args.audit_out),
    )


if __name__ == "__main__":
    main()
