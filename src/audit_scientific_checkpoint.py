from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .biological_replicates import (
    _resolve_column,
    attach_biological_replicates,
    canonical_replicate_column,
)
from .utils import ensure_dir, load_config


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sampled_file_fingerprint(path: Path) -> dict[str, Any]:
    stat = path.stat()
    sample_size = 1024 * 1024
    digest = hashlib.sha256()
    digest.update(str(stat.st_size).encode("ascii"))
    with path.open("rb") as handle:
        digest.update(handle.read(sample_size))
        if stat.st_size > sample_size:
            handle.seek(max(0, stat.st_size - sample_size))
            digest.update(handle.read(sample_size))
    return {
        "path": str(path.resolve()),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sampled_sha256": digest.hexdigest(),
        "sampled_sha256_scope": "file size plus first and last 1 MiB",
    }


def _resolved_config_fingerprint(config: dict[str, Any]) -> str:
    payload = json.dumps(
        config,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _implementation_fingerprint() -> dict[str, Any]:
    paths = sorted(Path("src").glob("*.py")) + sorted(Path("workflows").glob("Snakefile*"))
    digest = hashlib.sha256()
    included: list[str] = []
    for path in paths:
        relative = path.as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(path.read_bytes())
        included.append(relative)
    return {
        "sha256": digest.hexdigest(),
        "scope": "tracked Python modules and Snakemake workflow definitions",
        "n_files": len(included),
    }


def _check(
    checks: list[dict[str, Any]],
    *,
    name: str,
    passed: bool,
    detail: str,
    severity: str = "error",
) -> None:
    checks.append(
        {
            "name": name,
            "passed": bool(passed),
            "severity": severity,
            "detail": detail,
        }
    )


def _required_obs_columns(config: dict[str, Any]) -> list[str]:
    columns = {canonical_replicate_column(config)}
    unit_section = config.get("biological_units", {})
    for key in ("subject_col", "sample_unit_col", "technical_library_col"):
        value = unit_section.get(key)
        if value:
            columns.add(str(value))
    for section_name in ("composition_age", "signature_age", "age_prediction"):
        section = config.get(section_name, {})
        for key in ("age_col", "celltype_col"):
            value = section.get(key)
            if value:
                columns.add(str(value))
        columns.update(str(value) for value in section.get("covariate_cols", []))
    return sorted(columns)


def _stage_retention(config: dict[str, Any], checks: list[dict[str, Any]]) -> pd.DataFrame:
    stages = config.get("scientific_audit", {}).get("stage_h5ads", {})
    if not isinstance(stages, dict) or not stages:
        raise ValueError("scientific_audit.stage_h5ads must define ordered stages.")

    rows: list[dict[str, Any]] = []
    previous_count: int | None = None
    previous_stage: str | None = None
    for stage, value in stages.items():
        path = Path(str(value))
        if not path.is_file():
            raise FileNotFoundError(f"Configured audit stage is missing: {path}")
        stage_data = ad.read_h5ad(path, backed="r")
        try:
            count = int(stage_data.n_obs)
        finally:
            stage_data.file.close()
        input_count = count if previous_count is None else previous_count
        retained = float(count / input_count) if input_count else np.nan
        rows.append(
            {
                "scope": "workflow_stage",
                "stage": str(stage),
                "stratum": "all_cells",
                "level": "all",
                "input_cells": input_count,
                "retained_cells": count,
                "retention_fraction": retained,
                "decision": "configured workflow stage",
            }
        )
        if previous_count is not None:
            _check(
                checks,
                name=f"nonincreasing_cells__{stage}",
                passed=count <= previous_count,
                detail=(f"{previous_stage}={previous_count:,}; {stage}={count:,}"),
            )
        previous_count = count
        previous_stage = str(stage)
    return pd.DataFrame(rows)


def _retention_by_stratum(config: dict[str, Any], mapping_path: Path) -> pd.DataFrame:
    longitudinal = "biological_units" in config
    section = config.get("biological_units", config.get("biological_replicates", {}))
    source_path = Path(str(section["source_table_path"]))
    separator = str(section.get("source_table_sep", ","))
    header = list(pd.read_csv(source_path, sep=separator, nrows=0).columns)
    requests = {
        "cell_id": str(section.get("source_table_join_key", "cell_id")),
        "age": str(section.get("age_col", "age")),
        "sex": str(section.get("sex_col", "sex")),
        "batch": str(section.get("batch_col", "batch")),
    }
    if longitudinal:
        requests.update(
            {
                str(section.get("subject_col", "subject_id")): str(
                    section.get("subject_source_col", "donor_id")
                ),
                str(section.get("sample_unit_col", "sample_unit_id")): str(
                    section.get("sample_unit_source_col", "tube_id")
                ),
                str(section.get("technical_library_col", "technical_library_id")): str(
                    section.get("technical_library_source_col", "file_name")
                ),
            }
        )
    else:
        requests["biological_replicate_id"] = str(section.get("source_col", "tube_id"))
    resolved = {
        output: _resolve_column(header, requested) for output, requested in requests.items()
    }
    missing = [name for name, column in resolved.items() if column is None]
    if missing:
        raise KeyError(f"Source metadata cannot support retention strata: {missing}.")

    source = pd.read_csv(
        source_path,
        sep=separator,
        usecols=list(dict.fromkeys(resolved.values())),
        low_memory=False,
    ).rename(columns={value: key for key, value in resolved.items()})
    selected = pd.read_csv(mapping_path, low_memory=False)
    required_selected = set(requests)
    absent = sorted(required_selected.difference(selected.columns))
    if absent:
        raise KeyError(f"Replicate mapping lacks retention fields: {absent}.")

    for frame in (source, selected):
        frame["age"] = pd.to_numeric(frame["age"], errors="coerce")
        frame["age_decade"] = (np.floor(frame["age"] / 10.0) * 10.0).astype("Int64")

    rows: list[dict[str, Any]] = []
    strata = (
        [
            str(section.get("subject_col", "subject_id")),
            str(section.get("sample_unit_col", "sample_unit_id")),
            "sex",
            "batch",
            "age_decade",
        ]
        if longitudinal
        else ["biological_replicate_id", "sex", "batch", "age_decade"]
    )
    for stratum in strata:
        source_counts = source[stratum].astype("string").fillna("__MISSING__").value_counts()
        selected_counts = selected[stratum].astype("string").fillna("__MISSING__").value_counts()
        for level in sorted(set(source_counts.index).union(selected_counts.index)):
            input_cells = int(source_counts.get(level, 0))
            retained_cells = int(selected_counts.get(level, 0))
            rows.append(
                {
                    "scope": "source_to_checkpoint",
                    "stage": "annotated_checkpoint",
                    "stratum": stratum,
                    "level": str(level),
                    "input_cells": input_cells,
                    "retained_cells": retained_cells,
                    "retention_fraction": (
                        float(retained_cells / input_cells) if input_cells else np.nan
                    ),
                    "decision": "descriptive only; no automatic threshold",
                }
            )
    return pd.DataFrame(rows)


def _annotation_confidence(obs: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    section = config.get("scientific_audit", {})
    confidence_col = str(section.get("annotation_confidence_col", "cell_type_confidence"))
    celltype_col = str(config.get("composition_age", {}).get("celltype_col", "cell_type"))
    if confidence_col not in obs.columns or celltype_col not in obs.columns:
        return pd.DataFrame(
            columns=[
                "cell_type",
                "n_cells",
                "n_confidence_available",
                "confidence_q05",
                "confidence_median",
                "confidence_q95",
            ]
        )

    table = pd.DataFrame(
        {
            "cell_type": obs[celltype_col].astype("string"),
            "confidence": pd.to_numeric(obs[confidence_col], errors="coerce"),
        }
    )
    rows: list[dict[str, Any]] = []
    for cell_type, group in table.groupby("cell_type", observed=True):
        available = group["confidence"].dropna()
        rows.append(
            {
                "cell_type": str(cell_type),
                "n_cells": int(len(group)),
                "n_confidence_available": int(len(available)),
                "confidence_q05": float(available.quantile(0.05))
                if not available.empty
                else np.nan,
                "confidence_median": float(available.median()) if not available.empty else np.nan,
                "confidence_q95": float(available.quantile(0.95))
                if not available.empty
                else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("cell_type")


def _provenance_checks(
    adata: ad.AnnData,
    config: dict[str, Any],
    checks: list[dict[str, Any]],
) -> dict[str, Any]:
    required_fields = {
        "scvi_training_provenance": {
            "started_at_utc",
            "completed_at_utc",
            "input_path",
            "input_shape",
            "input_obs_names_sha256",
            "input_var_names_sha256",
            "config_path",
            "resolved_config_sha256",
            "seed",
            "layer",
            "categorical_covariates",
            "continuous_covariates",
            "allow_missing_covariates",
            "model_args",
            "training_args",
            "accelerator",
            "devices",
            "scvi_tools_version",
            "torch_version",
            "anndata_version",
            "scanpy_version",
        },
        "clustering_provenance": {
            "started_at_utc",
            "completed_at_utc",
            "input_path",
            "input_shape",
            "input_obs_names_sha256",
            "input_var_names_sha256",
            "config_path",
            "resolved_config_sha256",
            "seed",
            "n_jobs",
            "use_rep",
            "leiden_resolution",
            "neighbors_args",
            "umap_args",
            "leiden_args",
            "scanpy_version",
        },
        "celltypist_provenance": {
            "model_identifier",
            "resolved_model_path",
            "model_sha256",
            "model_source",
            "celltypist_version",
            "majority_voting",
            "chunk_size",
            "normalize_target_sum",
        },
    }
    observed: dict[str, Any] = {}
    for key, fields in required_fields.items():
        value = adata.uns.get(key)
        observed[key] = value
        missing_fields = (
            sorted(fields.difference(value.keys())) if isinstance(value, dict) else sorted(fields)
        )
        _check(
            checks,
            name=f"provenance__{key}",
            passed=isinstance(value, dict) and not missing_fields,
            detail=(
                "complete"
                if isinstance(value, dict) and not missing_fields
                else f"missing fields: {missing_fields}"
            ),
        )

    seed = int(config.get("run", {}).get("seed", 42))
    for key in ("scvi_training_provenance", "clustering_provenance"):
        value = observed.get(key)
        if isinstance(value, dict) and "seed" in value:
            _check(
                checks,
                name=f"seed_matches_config__{key}",
                passed=int(value["seed"]) == seed,
                detail=f"configured={seed}; observed={value['seed']}",
            )
    clustering = observed.get("clustering_provenance")
    clustering_config = config.get("clustering")
    if isinstance(clustering, dict) and isinstance(clustering_config, dict):
        from .cluster import _h5ad_safe_arguments

        expected = {
            "n_jobs": int(clustering_config["n_jobs"]),
            "use_rep": str(clustering_config["use_rep"]),
            "leiden_resolution": float(clustering_config["leiden_resolution"]),
            "neighbors_args": _h5ad_safe_arguments(clustering_config["neighbors_args"]),
            "umap_args": _h5ad_safe_arguments(clustering_config["umap_args"]),
            "leiden_args": _h5ad_safe_arguments(clustering_config["leiden_args"]),
        }
        mismatches = {
            key: {"configured": value, "observed": clustering.get(key)}
            for key, value in expected.items()
            if clustering.get(key) != value
        }
        _check(
            checks,
            name="clustering_parameters_match_config",
            passed=not mismatches,
            detail="matched" if not mismatches else f"mismatches={mismatches}",
        )
    return observed


def audit(
    *,
    config_path: Path,
    checkpoint_path: Path,
    mapping_path: Path,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    config = load_config(config_path)
    unit_section = config.get("biological_units", config.get("biological_replicates", {}))
    replicate_source_path = Path(str(unit_section["source_table_path"]))
    checks: list[dict[str, Any]] = []
    retention = _stage_retention(config, checks)
    retention = pd.concat(
        [retention, _retention_by_stratum(config, mapping_path)],
        ignore_index=True,
    )

    adata = ad.read_h5ad(checkpoint_path, backed="r")
    try:
        replicate_report = attach_biological_replicates(adata, config)
        obs = adata.obs.copy()
        required_obs = _required_obs_columns(config)
        missing_obs = sorted(set(required_obs).difference(obs.columns))
        _check(
            checks,
            name="required_obs_columns",
            passed=not missing_obs,
            detail=f"required={required_obs}; missing={missing_obs}",
        )
        missing_values = {
            column: int(obs[column].isna().sum())
            for column in required_obs
            if column in obs.columns and obs[column].isna().any()
        }
        _check(
            checks,
            name="required_obs_complete",
            passed=not missing_values,
            detail=f"missing values={missing_values}",
        )
        blank_values = {
            column: int(obs[column].astype("string").str.strip().eq("").sum())
            for column in required_obs
            if column in obs.columns and obs[column].astype("string").str.strip().eq("").any()
        }
        _check(
            checks,
            name="required_obs_nonblank",
            passed=not blank_values,
            detail=f"blank values={blank_values}",
        )
        age_columns = sorted(
            {
                str(config[section]["age_col"])
                for section in (
                    "composition_age",
                    "signature_age",
                    "age_prediction",
                )
                if config.get(section, {}).get("age_col")
            }
        )
        invalid_ages = {
            column: int(pd.to_numeric(obs[column], errors="coerce").isna().sum())
            for column in age_columns
            if column in obs.columns and pd.to_numeric(obs[column], errors="coerce").isna().any()
        }
        _check(
            checks,
            name="age_values_numeric",
            passed=not invalid_ages,
            detail=f"non-numeric or missing ages={invalid_ages}",
        )
        _check(
            checks,
            name="unique_cell_identifiers",
            passed=bool(adata.obs_names.is_unique),
            detail=f"n_obs={adata.n_obs:,}",
        )
        _check(
            checks,
            name="unique_gene_identifiers",
            passed=bool(adata.var_names.is_unique),
            detail=f"n_vars={adata.n_vars:,}",
        )
        latent_key = str(config.get("age_prediction", {}).get("latent_key", "X_scVI"))
        _check(
            checks,
            name="configured_latent_representation",
            passed=latent_key in adata.obsm,
            detail=f"required={latent_key}; available={sorted(adata.obsm.keys())}",
        )
        provenance = _provenance_checks(adata, config, checks)
        confidence = _annotation_confidence(obs, config)
        _check(
            checks,
            name="annotation_confidence_available",
            passed=not confidence.empty,
            detail=(
                f"summarized {len(confidence)} cell types"
                if not confidence.empty
                else "cell-level annotation confidence is unavailable"
            ),
            severity="warning",
        )
        matrix_state = {
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "x_type": type(adata.X).__name__,
            "x_dtype": str(adata.X.dtype),
            "validation_scope": (
                "structural only; pseudobulk aggregation performs the exact "
                "count-state and conservation checks"
            ),
        }
    finally:
        adata.file.close()

    _check(
        checks,
        name="placeholder_outputs_disabled",
        passed=not bool(config.get("run", {}).get("allow_placeholder_outputs", False)),
        detail=("maintained study profiles must fail on analytical errors"),
    )
    errors = [item for item in checks if item["severity"] == "error" and not item["passed"]]
    report = {
        "passed": not errors,
        "scope": (
            "scientific checkpoint compatibility and descriptive QC; not biological validation"
        ),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "config": {
            "path": str(config_path.resolve()),
            "file_sha256": _sha256_file(config_path),
            "resolved_sha256": _resolved_config_fingerprint(config),
        },
        "checkpoint": _sampled_file_fingerprint(checkpoint_path),
        "replicate_source_metadata": _sampled_file_fingerprint(replicate_source_path),
        "replicate_mapping": _sampled_file_fingerprint(mapping_path),
        "implementation": _implementation_fingerprint(),
        "replicate_attachment": replicate_report,
        "matrix_state": matrix_state,
        "qc_decisions": {
            "qc_enabled": bool(config.get("qc", {}).get("enabled", True)),
            "doublet_removal_enabled": bool(config.get("doublets", {}).get("enabled", True)),
            "interpretation": (
                "configuration record only; disabled stages require human "
                "review and are not treated as evidence of clean data"
            ),
        },
        "provenance": provenance,
        "checks": checks,
        "failed_checks": [item["name"] for item in errors],
    }
    return report, retention, confidence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--retention-out", required=True)
    parser.add_argument("--confidence-out", required=True)
    args = parser.parse_args()

    report, retention, confidence = audit(
        config_path=Path(args.config),
        checkpoint_path=Path(args.inp),
        mapping_path=Path(args.mapping),
    )
    report_path = Path(args.report)
    retention_path = Path(args.retention_out)
    confidence_path = Path(args.confidence_out)
    for parent in (
        report_path.parent,
        retention_path.parent,
        confidence_path.parent,
    ):
        ensure_dir(parent)
    retention.to_csv(retention_path, index=False)
    confidence.to_csv(confidence_path, index=False)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    if not report["passed"]:
        raise RuntimeError(
            "Scientific checkpoint audit failed: " + ", ".join(report["failed_checks"])
        )
    print(
        f"[audit_scientific_checkpoint] passed {len(report['checks'])} checks for {args.inp}",
        flush=True,
    )


if __name__ == "__main__":
    main()
