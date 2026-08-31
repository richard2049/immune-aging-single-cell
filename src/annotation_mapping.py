from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pandas as pd

MAPPING_COLUMNS = (
    "raw_label",
    "analysis_label",
    "disposition",
    "rationale",
    "review_status",
    "mapping_version",
)
REQUIRED_MAPPING_COLUMNS = set(MAPPING_COLUMNS)
VALID_DISPOSITIONS = {"primary", "exploratory", "unresolved"}
VALID_REVIEW_STATUSES = {"pending_review", "approved", "rejected"}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def annotation_config(cfg: dict[str, Any]) -> dict[str, Any]:
    section = cfg.get("annotation_qualification", {})
    if not isinstance(section, dict):
        raise TypeError("annotation_qualification must be a mapping.")
    return section


def load_label_mapping(
    path: Path,
    *,
    observed_labels: set[str] | None = None,
) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Annotation label mapping not found: {path}")
    mapping = pd.read_csv(path, dtype="string", keep_default_na=False)
    missing_columns = sorted(REQUIRED_MAPPING_COLUMNS.difference(mapping.columns))
    if missing_columns:
        raise KeyError("Annotation mapping is missing columns: " + ", ".join(missing_columns))

    mapping = mapping[list(MAPPING_COLUMNS)].copy()
    for column in MAPPING_COLUMNS:
        mapping[column] = mapping[column].astype("string").str.strip()
    if mapping.empty:
        raise ValueError("Annotation mapping is empty.")
    if mapping["raw_label"].eq("").any():
        raise ValueError("Annotation mapping contains an empty raw_label.")
    if mapping["raw_label"].duplicated().any():
        duplicates = sorted(
            mapping.loc[mapping["raw_label"].duplicated(False), "raw_label"].unique()
        )
        raise ValueError(
            "Annotation mapping contains duplicate raw labels: " + ", ".join(duplicates)
        )
    invalid_dispositions = sorted(set(mapping["disposition"]).difference(VALID_DISPOSITIONS))
    if invalid_dispositions:
        raise ValueError("Invalid annotation dispositions: " + ", ".join(invalid_dispositions))
    invalid_statuses = sorted(set(mapping["review_status"]).difference(VALID_REVIEW_STATUSES))
    if invalid_statuses:
        raise ValueError("Invalid annotation review statuses: " + ", ".join(invalid_statuses))
    if mapping["rationale"].eq("").any() or mapping["mapping_version"].eq("").any():
        raise ValueError("Every annotation mapping row requires rationale and mapping_version.")
    if mapping.loc[mapping["disposition"].ne("unresolved"), "analysis_label"].eq("").any():
        raise ValueError("Primary and exploratory rows require an analysis_label.")

    versions = mapping["mapping_version"].unique()
    if len(versions) != 1:
        raise ValueError("Annotation mapping must contain exactly one mapping_version.")
    if observed_labels is not None:
        mapped_labels = set(mapping["raw_label"].astype(str))
        missing = sorted(observed_labels.difference(mapped_labels))
        extra = sorted(mapped_labels.difference(observed_labels))
        if missing or extra:
            raise ValueError(
                "Annotation mapping does not match observed labels; "
                f"missing={missing}; extra={extra}"
            )
    return mapping.sort_values("raw_label").reset_index(drop=True)


def mapping_is_approved(mapping: pd.DataFrame) -> bool:
    return bool(mapping["review_status"].eq("approved").all())


def attach_analysis_cell_types(adata, cfg: dict[str, Any]) -> dict[str, Any]:
    section = annotation_config(cfg)
    if not bool(section.get("enabled", False)):
        return {"enabled": False}

    raw_col = str(section.get("raw_label_col", "cell_type"))
    analysis_col = str(section.get("analysis_label_col", "cell_type_analysis"))
    disposition_col = str(section.get("disposition_col", "cell_type_analysis_disposition"))
    if raw_col not in adata.obs.columns:
        raise KeyError(f"Raw annotation column {raw_col!r} is absent from AnnData.obs.")

    mapping_path = Path(str(section.get("mapping_path", "")))
    observed_labels = set(adata.obs[raw_col].dropna().astype(str).unique())
    mapping = load_label_mapping(mapping_path, observed_labels=observed_labels)
    if bool(section.get("require_approved_mapping", True)) and not mapping_is_approved(mapping):
        pending = int(mapping["review_status"].ne("approved").sum())
        raise RuntimeError(
            f"Annotation mapping is not approved: {pending} rows remain pending or rejected."
        )

    lookup = mapping.set_index("raw_label")
    raw = adata.obs[raw_col].astype("string")
    disposition = raw.map(lookup["disposition"])
    analysis = raw.map(lookup["analysis_label"])
    analysis = analysis.where(disposition.eq("primary"), pd.NA)
    adata.obs[analysis_col] = analysis.to_numpy()
    adata.obs[disposition_col] = disposition.to_numpy()
    return {
        "enabled": True,
        "mapping_path": str(mapping_path),
        "mapping_sha256": file_sha256(mapping_path),
        "mapping_version": str(mapping["mapping_version"].iloc[0]),
        "raw_label_col": raw_col,
        "analysis_label_col": analysis_col,
        "primary_cells": int(analysis.notna().sum()),
        "nonprimary_cells": int(analysis.isna().sum()),
        "primary_analysis_labels": int(analysis.nunique(dropna=True)),
    }
