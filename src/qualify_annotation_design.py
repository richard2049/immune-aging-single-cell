from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .annotation_mapping import (
    annotation_config,
    file_sha256,
    load_label_mapping,
    mapping_is_approved,
)
from .scvi_train import _resolved_config_sha256
from .utils import ensure_dir, load_config
from .validate_scvi_checkpoint import _sampled_file_fingerprint


def _atomic_csv(table: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    table.to_csv(temporary, index=False)
    os.replace(temporary, path)


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def _align_units(obs_names: pd.Index, path: Path) -> pd.DataFrame:
    required = [
        "cell_id",
        "subject_id",
        "sample_unit_id",
        "technical_library_id",
        "age",
        "sex",
        "batch",
    ]
    mapping = pd.read_csv(path, usecols=required, low_memory=False)
    if mapping["cell_id"].duplicated().any():
        raise ValueError("Cell-to-unit mapping contains duplicate cell IDs.")
    observed = obs_names.astype(str).to_numpy()
    mapped = mapping["cell_id"].astype(str).to_numpy()
    if not np.array_equal(observed, mapped):
        mapping = (
            mapping.assign(_cell_id_key=mapping["cell_id"].astype(str))
            .set_index("_cell_id_key")
            .reindex(observed)
        )
        if mapping["sample_unit_id"].isna().any():
            raise ValueError("Cell-to-unit mapping does not cover the checkpoint exactly.")
        mapping = mapping.reset_index(drop=True)
        mapping["cell_id"] = observed
        mapping = mapping[required]
    return mapping


def _review_markdown(
    *,
    label_mapping: pd.DataFrame,
    evidence: pd.DataFrame,
    design: pd.DataFrame,
    n_cells: int,
) -> str:
    proposed = evidence.loc[evidence["label_level"].eq("proposed_analysis")].copy()
    proposed = proposed.sort_values("n_cells", ascending=False)
    proposed_design = design.loc[design["label_level"].eq("proposed_analysis")]
    passing = (
        proposed_design.groupby("method", observed=True)["passes_current_contract"]
        .sum()
        .astype(int)
        .to_dict()
    )
    raw_label_count = int(evidence.loc[evidence["label_level"].eq("raw"), "label"].nunique())
    status_counts = ", ".join(
        f"{status}={count}"
        for status, count in label_mapping["review_status"].value_counts().items()
    )
    approval = "approved" if mapping_is_approved(label_mapping) else "blocked pending human review"
    if mapping_is_approved(label_mapping):
        review_instruction = [
            "The approved mapping reflects the completed human review. Retain the mapping",
            "rationale, cluster cross-tabulation, and descriptive marker table with this",
            "qualification record for provenance.",
        ]
    else:
        review_instruction = [
            "Review the mapping rationale, cluster cross-tabulation, and descriptive marker",
            "table before changing all mapping rows to `approved`.",
        ]
    lines = [
        "# Annotation And Design Qualification Review",
        "",
        "This report is outcome-blind. It contains annotation support and design",
        "eligibility only; no age associations, effect directions, or significance",
        "results were inspected or generated.",
        "",
        "## Review Status",
        "",
        f"- Cells assessed: {n_cells:,}",
        f"- Raw labels: {raw_label_count}",
        f"- Proposed primary labels: {len(proposed)}",
        f"- Mapping rows by review status: {status_counts}",
        f"- Downstream approval: {approval}",
        "",
        "## Proposed Primary Populations",
        "",
        "| Analysis label | Cells | Sample units | Subjects | Batches | Median confidence |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in proposed.itertuples(index=False):
        lines.append(
            f"| {row.label} | {row.n_cells:,} | {row.n_sample_units_detected} | "
            f"{row.n_subjects_detected} | {row.n_batches_detected} | "
            f"{row.confidence_median:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Method-Specific Eligibility",
            "",
            *[
                f"- `{method}`: {count} proposed populations pass the current structural contract."
                for method, count in sorted(passing.items())
            ],
            "",
            "Passing these checks establishes technical support and estimability under",
            "the configured design. It does not validate biological identity or approve",
            "an immune-ageing claim.",
            *review_instruction,
            "",
        ]
    )
    return "\n".join(lines)


def _design_diagnostics(metadata: pd.DataFrame) -> dict[str, Any]:
    if metadata.empty:
        return {
            "design_columns": 0,
            "design_rank": 0,
            "residual_df": 0,
            "condition_number": np.nan,
            "design_full_rank": False,
        }
    age = pd.to_numeric(metadata["age"], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(age).all():
        raise ValueError("Design metadata contain missing or non-numeric age values.")
    age_sd = float(np.std(age, ddof=0))
    age_scaled = (age - float(np.mean(age))) / age_sd if age_sd > 0 else np.zeros(len(age))
    nuisance = pd.get_dummies(
        metadata[["sex", "batch"]].astype("string"),
        drop_first=True,
        dtype=float,
    )
    matrix = np.column_stack([np.ones(len(metadata)), age_scaled, nuisance.to_numpy()])
    rank = int(np.linalg.matrix_rank(matrix))
    return {
        "design_columns": int(matrix.shape[1]),
        "design_rank": rank,
        "residual_df": int(matrix.shape[0] - rank),
        "condition_number": float(np.linalg.cond(matrix)),
        "design_full_rank": bool(rank == matrix.shape[1]),
    }


def _label_evidence(
    obs: pd.DataFrame,
    *,
    label_col: str,
    label_level: str,
    low_confidence_reference: float,
) -> pd.DataFrame:
    valid = obs.dropna(subset=[label_col]).copy()
    valid[label_col] = valid[label_col].astype(str)
    rows: list[dict[str, Any]] = []
    for label, group in valid.groupby(label_col, observed=True, sort=True):
        cluster_counts = group["leiden"].astype(str).value_counts()
        batch_counts = group["batch"].astype(str).value_counts()
        confidence = pd.to_numeric(group["cell_type_confidence"], errors="coerce")
        rows.append(
            {
                "label_level": label_level,
                "label": str(label),
                "n_cells": int(len(group)),
                "n_sample_units_detected": int(group["sample_unit_id"].nunique()),
                "n_subjects_detected": int(group["subject_id"].nunique()),
                "n_batches_detected": int(group["batch"].nunique()),
                "age_min": float(group["age"].min()),
                "age_max": float(group["age"].max()),
                "age_span": float(group["age"].max() - group["age"].min()),
                "confidence_q05": float(confidence.quantile(0.05)),
                "confidence_median": float(confidence.median()),
                "confidence_q95": float(confidence.quantile(0.95)),
                "low_confidence_fraction": float(confidence.lt(low_confidence_reference).mean()),
                "dominant_leiden": str(cluster_counts.index[0]),
                "dominant_leiden_fraction": float(cluster_counts.iloc[0] / len(group)),
                "n_leiden_clusters": int(len(cluster_counts)),
                "dominant_batch": str(batch_counts.index[0]),
                "dominant_batch_fraction": float(batch_counts.iloc[0] / len(group)),
            }
        )
    return pd.DataFrame(rows)


def _method_contracts(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    composition = cfg["composition_age"]
    signature = cfg["signature_age"]
    prediction = cfg["age_prediction"]
    pseudobulk = cfg["pseudobulk_de"]
    return [
        {
            "method": "composition",
            "min_cells_per_profile": 1,
            "min_sample_units": 0,
            "min_subjects": int(composition["min_donors_per_celltype"]),
            "min_age_span": 0.0,
            "min_residual_df": 1,
            "requires_full_rank": True,
            "structural_zeros": True,
        },
        {
            "method": "signature",
            "min_cells_per_profile": int(signature["min_cells_per_group"]),
            "min_sample_units": int(signature["min_donors_per_celltype"]),
            "min_subjects": int(signature["min_donors_per_celltype"]),
            "min_age_span": float(signature["min_age_span"]),
            "min_residual_df": 1,
            "requires_full_rank": True,
            "structural_zeros": False,
        },
        {
            "method": "age_prediction_figure",
            "min_cells_per_profile": int(prediction["min_cells_per_group"]),
            "min_sample_units": int(prediction["min_groups_for_mae_plot"]),
            "min_subjects": int(prediction["cv_folds"]),
            "min_age_span": 0.0,
            "min_residual_df": 0,
            "requires_full_rank": False,
            "structural_zeros": False,
        },
        {
            "method": "pseudobulk",
            "min_cells_per_profile": int(pseudobulk["min_cells_per_pseudobulk"]),
            "min_sample_units": int(pseudobulk["min_replicates_per_celltype"]),
            "min_subjects": int(pseudobulk["min_subjects_per_celltype"]),
            "min_age_span": float(pseudobulk["min_age_span_years"]),
            "min_residual_df": int(pseudobulk["min_residual_df"]),
            "requires_full_rank": True,
            "structural_zeros": False,
        },
    ]


def _design_eligibility(
    obs: pd.DataFrame,
    unit_metadata: pd.DataFrame,
    *,
    label_col: str,
    label_level: str,
    cfg: dict[str, Any],
) -> pd.DataFrame:
    valid = obs.dropna(subset=[label_col]).copy()
    valid[label_col] = valid[label_col].astype(str)
    profiles = (
        valid.groupby([label_col, "sample_unit_id"], observed=True, sort=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .merge(unit_metadata, on="sample_unit_id", validate="many_to_one")
    )
    rows: list[dict[str, Any]] = []
    for contract in _method_contracts(cfg):
        for label, group in profiles.groupby(label_col, observed=True, sort=True):
            detected_subjects = int(group["subject_id"].nunique())
            if contract["structural_zeros"]:
                eligible = unit_metadata.copy()
                n_profiles = int(len(unit_metadata))
                n_subjects = int(unit_metadata["subject_id"].nunique())
            else:
                eligible = group.loc[
                    group["n_cells"].ge(contract["min_cells_per_profile"]),
                    ["sample_unit_id", "subject_id", "age", "sex", "batch"],
                ].copy()
                n_profiles = int(eligible["sample_unit_id"].nunique())
                n_subjects = int(eligible["subject_id"].nunique())
            age_span = (
                float(eligible["age"].max() - eligible["age"].min()) if not eligible.empty else 0.0
            )
            design = _design_diagnostics(eligible)
            support_subjects = detected_subjects if contract["structural_zeros"] else n_subjects
            passed = (
                n_profiles >= contract["min_sample_units"]
                and support_subjects >= contract["min_subjects"]
                and age_span >= contract["min_age_span"]
                and design["residual_df"] >= contract["min_residual_df"]
                and (design["design_full_rank"] if contract["requires_full_rank"] else True)
            )
            rows.append(
                {
                    "label_level": label_level,
                    "label": str(label),
                    "method": contract["method"],
                    "eligible_sample_units": n_profiles,
                    "eligible_subjects": n_subjects,
                    "detected_subjects": detected_subjects,
                    "eligible_batches": int(eligible["batch"].nunique()),
                    "eligible_age_span": age_span,
                    "min_cells_per_profile": contract["min_cells_per_profile"],
                    "required_sample_units": contract["min_sample_units"],
                    "required_subjects": contract["min_subjects"],
                    "required_age_span": contract["min_age_span"],
                    "required_residual_df": contract["min_residual_df"],
                    **design,
                    "passes_current_contract": bool(passed),
                }
            )
    return pd.DataFrame(rows)


def _cluster_crosstab(obs: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for level, column in (("raw", "cell_type"), ("proposed_analysis", "cell_type_analysis")):
        valid = obs.dropna(subset=[column])
        table = (
            valid.groupby(["leiden", column], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .rename(columns={column: "label"})
        )
        table.insert(0, "label_level", level)
        table["cluster_fraction"] = table["n_cells"] / table.groupby("leiden")["n_cells"].transform(
            "sum"
        )
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def _marker_summary_tables(
    adata,
    label_sets: dict[str, pd.Series],
    panel_by_label: dict[str, dict[str, str]],
    marker_panels: dict[str, list[str]],
    *,
    chunk_size: int,
) -> dict[str, pd.DataFrame]:
    requested = sorted({gene for genes in marker_panels.values() for gene in genes})
    if not requested:
        raise ValueError("Annotation qualification requires at least one review marker.")
    var_lookup = {str(name).upper(): index for index, name in enumerate(adata.var_names)}
    present = [gene for gene in requested if gene.upper() in var_lookup]
    present_indices = [var_lookup[gene.upper()] for gene in present]
    present_position = {gene: index for index, gene in enumerate(present)}
    states: dict[str, dict[str, Any]] = {}
    for level, labels in label_sets.items():
        groups = sorted(labels.dropna().astype(str).unique())
        states[level] = {
            "labels": labels,
            "group_index": {group: index for index, group in enumerate(groups)},
            "sums": np.zeros((len(groups), len(present)), dtype=np.float64),
            "detected": np.zeros((len(groups), len(present)), dtype=np.int64),
            "n_cells": np.zeros(len(groups), dtype=np.int64),
        }

    total_chunks = int(np.ceil(adata.n_obs / chunk_size))
    for chunk_number, start in enumerate(range(0, adata.n_obs, chunk_size), start=1):
        end = min(start + chunk_size, adata.n_obs)
        if not present_indices:
            continue
        matrix = adata.X[start:end]
        matrix = sparse.csr_matrix(matrix) if not sparse.issparse(matrix) else matrix.tocsr()
        matrix = matrix[:, present_indices]
        for state in states.values():
            labels = state["labels"].iloc[start:end].astype("string")
            for group in labels.dropna().unique():
                mask = labels.eq(group).fillna(False).to_numpy(dtype=bool)
                selected = matrix[mask]
                index = state["group_index"][str(group)]
                state["n_cells"][index] += int(mask.sum())
                state["sums"][index] += np.asarray(selected.sum(axis=0)).ravel()
                state["detected"][index] += np.asarray(selected.getnnz(axis=0)).ravel()
        if chunk_number == 1 or chunk_number % 10 == 0 or chunk_number == total_chunks:
            print(
                f"[qualify_annotation_design] marker scan chunk {chunk_number}/{total_chunks}",
                flush=True,
            )

    outputs: dict[str, pd.DataFrame] = {}
    for level, state in states.items():
        rows: list[dict[str, Any]] = []
        for group, group_position in state["group_index"].items():
            panel_label = panel_by_label.get(level, {}).get(group, "")
            expected = set(marker_panels.get(panel_label, []))
            denominator = int(state["n_cells"][group_position])
            for gene in requested:
                if gene in present:
                    position = present_position[gene]
                    detection_fraction = (
                        float(state["detected"][group_position, position] / denominator)
                        if denominator
                        else np.nan
                    )
                    mean_raw_count = (
                        float(state["sums"][group_position, position] / denominator)
                        if denominator
                        else np.nan
                    )
                else:
                    detection_fraction = np.nan
                    mean_raw_count = np.nan
                rows.append(
                    {
                        "label_level": level,
                        "label": group,
                        "marker_panel_label": panel_label,
                        "gene": gene,
                        "expected_marker": gene in expected,
                        "gene_present": gene in present,
                        "n_cells": denominator,
                        "detection_fraction": detection_fraction,
                        "mean_raw_count": mean_raw_count,
                        "expression_state": "raw counts; descriptive marker evidence only",
                    }
                )
        outputs[level] = pd.DataFrame(rows)
    return outputs


def _mapping_review_table(
    label_mapping: pd.DataFrame,
    raw_evidence: pd.DataFrame,
    raw_design: pd.DataFrame,
    raw_marker: pd.DataFrame,
) -> pd.DataFrame:
    review = label_mapping.merge(
        raw_evidence.drop(columns="label_level").rename(columns={"label": "raw_label"}),
        on="raw_label",
        how="left",
        validate="one_to_one",
    )
    if review["n_cells"].isna().any():
        raise RuntimeError("Raw-label evidence does not cover the annotation mapping.")

    for method, group in raw_design.groupby("method", observed=True):
        fields = group[
            [
                "label",
                "eligible_sample_units",
                "eligible_subjects",
                "eligible_batches",
                "eligible_age_span",
                "design_full_rank",
                "residual_df",
                "passes_current_contract",
            ]
        ].rename(
            columns={
                "label": "raw_label",
                **{
                    column: f"{method}_{column}"
                    for column in (
                        "eligible_sample_units",
                        "eligible_subjects",
                        "eligible_batches",
                        "eligible_age_span",
                        "design_full_rank",
                        "residual_df",
                        "passes_current_contract",
                    )
                },
            }
        )
        review = review.merge(fields, on="raw_label", how="left", validate="one_to_one")

    marker_rows: list[dict[str, Any]] = []
    for raw_label, group in raw_marker.groupby("raw_label", observed=True, sort=True):
        expected = group.loc[group["expected_marker"]]
        detection = "; ".join(
            f"{row.gene}={row.detection_fraction:.3f}"
            for row in expected.loc[expected["gene_present"]].itertuples(index=False)
        )
        marker_rows.append(
            {
                "raw_label": raw_label,
                "marker_panel_label": str(group["marker_panel_label"].iloc[0]),
                "expected_marker_count": int(len(expected)),
                "expected_markers_present": int(expected["gene_present"].sum()),
                "expected_marker_detection": detection,
            }
        )
    review = review.merge(
        pd.DataFrame(marker_rows),
        on="raw_label",
        how="left",
        validate="one_to_one",
    )
    return review.sort_values(["disposition", "analysis_label", "raw_label"]).reset_index(drop=True)


def qualify(
    *,
    config_path: Path,
    checkpoint_path: Path,
    unit_mapping_path: Path,
    label_mapping_path: Path,
    evidence_path: Path,
    design_path: Path,
    cluster_path: Path,
    marker_path: Path,
    raw_marker_path: Path,
    review_table_path: Path,
    review_path: Path,
    report_path: Path,
) -> dict[str, Any]:
    cfg = load_config(config_path)
    section = annotation_config(cfg)
    raw_col = str(section.get("raw_label_col", "cell_type"))
    confidence_col = str(section.get("confidence_col", "cell_type_confidence"))
    cluster_col = str(section.get("cluster_col", "leiden"))
    analysis_col = str(section.get("analysis_label_col", "cell_type_analysis"))
    low_confidence_reference = float(section.get("low_confidence_reference", 0.5))

    adata = ad.read_h5ad(checkpoint_path, backed="r")
    try:
        required_obs = [raw_col, confidence_col, cluster_col]
        missing = [column for column in required_obs if column not in adata.obs.columns]
        if missing:
            raise KeyError("Checkpoint is missing annotation columns: " + ", ".join(missing))
        observed_labels = set(adata.obs[raw_col].dropna().astype(str).unique())
        label_mapping = load_label_mapping(
            label_mapping_path,
            observed_labels=observed_labels,
        )
        units = _align_units(adata.obs_names, unit_mapping_path)
        obs = pd.DataFrame(
            {
                "cell_type": adata.obs[raw_col].astype(str).to_numpy(),
                "cell_type_confidence": pd.to_numeric(
                    adata.obs[confidence_col], errors="coerce"
                ).to_numpy(),
                "leiden": adata.obs[cluster_col].astype(str).to_numpy(),
                "subject_id": units["subject_id"].astype(str).to_numpy(),
                "sample_unit_id": units["sample_unit_id"].astype(str).to_numpy(),
                "technical_library_id": units["technical_library_id"].astype(str).to_numpy(),
                "age": pd.to_numeric(units["age"], errors="coerce").to_numpy(),
                "sex": units["sex"].astype(str).to_numpy(),
                "batch": units["batch"].astype(str).to_numpy(),
            },
            index=adata.obs_names,
        )
        if obs.isna().any().any():
            raise ValueError("Annotation qualification inputs contain missing required values.")

        lookup = label_mapping.set_index("raw_label")
        disposition = obs["cell_type"].map(lookup["disposition"])
        proposed = obs["cell_type"].map(lookup["analysis_label"])
        obs[analysis_col] = proposed.where(disposition.eq("primary"), pd.NA)
        obs["mapping_disposition"] = disposition

        raw_evidence = _label_evidence(
            obs,
            label_col="cell_type",
            label_level="raw",
            low_confidence_reference=low_confidence_reference,
        )
        proposed_evidence = _label_evidence(
            obs,
            label_col=analysis_col,
            label_level="proposed_analysis",
            low_confidence_reference=low_confidence_reference,
        )
        evidence = pd.concat([raw_evidence, proposed_evidence], ignore_index=True)

        unit_metadata = units[
            ["sample_unit_id", "subject_id", "age", "sex", "batch"]
        ].drop_duplicates("sample_unit_id")
        raw_design = _design_eligibility(
            obs,
            unit_metadata,
            label_col="cell_type",
            label_level="raw",
            cfg=cfg,
        )
        proposed_design = _design_eligibility(
            obs,
            unit_metadata,
            label_col=analysis_col,
            label_level="proposed_analysis",
            cfg=cfg,
        )
        design = pd.concat([raw_design, proposed_design], ignore_index=True)
        cluster = _cluster_crosstab(obs)
        marker_panels = {
            str(label): [str(gene) for gene in genes]
            for label, genes in section.get("marker_panels", {}).items()
        }
        raw_panel_lookup = {
            str(row.raw_label): str(row.analysis_label)
            if str(row.analysis_label) in marker_panels
            else ""
            for row in label_mapping.itertuples(index=False)
        }
        marker_tables = _marker_summary_tables(
            adata,
            {
                "raw": obs["cell_type"],
                "proposed_analysis": obs[analysis_col],
            },
            {
                "raw": raw_panel_lookup,
                "proposed_analysis": {
                    str(label): str(label)
                    for label in obs[analysis_col].dropna().astype(str).unique()
                },
            },
            marker_panels,
            chunk_size=int(section.get("marker_chunk_size", 10_000)),
        )
        marker = marker_tables["proposed_analysis"].rename(columns={"label": "analysis_label"})
        raw_marker = marker_tables["raw"].rename(columns={"label": "raw_label"})
        raw_marker = raw_marker.merge(
            label_mapping[["raw_label", "analysis_label", "disposition"]],
            on="raw_label",
            how="left",
            validate="many_to_one",
        )
        review_table = _mapping_review_table(
            label_mapping,
            raw_evidence,
            raw_design,
            raw_marker,
        )
    finally:
        adata.file.close()

    _atomic_csv(evidence, evidence_path)
    _atomic_csv(design, design_path)
    _atomic_csv(cluster, cluster_path)
    _atomic_csv(marker, marker_path)
    _atomic_csv(raw_marker, raw_marker_path)
    _atomic_csv(review_table, review_table_path)
    ensure_dir(review_path.parent)
    review_temporary = review_path.with_name(review_path.name + ".tmp")
    review_temporary.write_text(
        _review_markdown(
            label_mapping=label_mapping,
            evidence=evidence,
            design=design,
            n_cells=len(obs),
        ),
        encoding="utf-8",
    )
    os.replace(review_temporary, review_path)
    review_counts = {
        str(key): int(value) for key, value in label_mapping["review_status"].value_counts().items()
    }
    disposition_counts = {
        str(key): int(value) for key, value in label_mapping["disposition"].value_counts().items()
    }
    approved = mapping_is_approved(label_mapping)
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "outcome-blind annotation and design qualification; no age associations",
        "status": "approved" if approved else "pending_review",
        "structural_checks_passed": True,
        "approved_for_downstream": approved,
        "config": {
            "path": str(config_path.resolve()),
            "resolved_sha256": _resolved_config_sha256(cfg),
        },
        "checkpoint": _sampled_file_fingerprint(checkpoint_path),
        "unit_mapping": _sampled_file_fingerprint(unit_mapping_path),
        "label_mapping": {
            "path": str(label_mapping_path.resolve()),
            "sha256": file_sha256(label_mapping_path),
            "version": str(label_mapping["mapping_version"].iloc[0]),
            "review_status_counts": review_counts,
            "disposition_counts": disposition_counts,
        },
        "summary": {
            "n_cells": int(len(obs)),
            "raw_labels": int(obs["cell_type"].nunique()),
            "proposed_primary_analysis_labels": int(obs[analysis_col].nunique(dropna=True)),
            "proposed_primary_cells": int(obs[analysis_col].notna().sum()),
            "nonprimary_cells": int(obs[analysis_col].isna().sum()),
            "passing_design_rows": {
                method: int(
                    design.loc[
                        design["label_level"].eq("proposed_analysis") & design["method"].eq(method),
                        "passes_current_contract",
                    ].sum()
                )
                for method in sorted(design["method"].unique())
            },
        },
        "outputs": {
            "evidence": _sampled_file_fingerprint(evidence_path),
            "design": _sampled_file_fingerprint(design_path),
            "cluster": _sampled_file_fingerprint(cluster_path),
            "marker": _sampled_file_fingerprint(marker_path),
            "raw_marker": _sampled_file_fingerprint(raw_marker_path),
            "mapping_review_table": _sampled_file_fingerprint(review_table_path),
            "review": _sampled_file_fingerprint(review_path),
        },
    }
    _atomic_json(report, report_path)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--inp", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--label-map", required=True)
    parser.add_argument("--evidence-out", required=True)
    parser.add_argument("--design-out", required=True)
    parser.add_argument("--cluster-out", required=True)
    parser.add_argument("--marker-out", required=True)
    parser.add_argument("--raw-marker-out", required=True)
    parser.add_argument("--review-table-out", required=True)
    parser.add_argument("--review-out", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    report = qualify(
        config_path=Path(args.config),
        checkpoint_path=Path(args.inp),
        unit_mapping_path=Path(args.mapping),
        label_mapping_path=Path(args.label_map),
        evidence_path=Path(args.evidence_out),
        design_path=Path(args.design_out),
        cluster_path=Path(args.cluster_out),
        marker_path=Path(args.marker_out),
        raw_marker_path=Path(args.raw_marker_out),
        review_table_path=Path(args.review_table_out),
        review_path=Path(args.review_out),
        report_path=Path(args.report),
    )
    print(
        "[qualify_annotation_design] "
        f"status={report['status']}; raw_labels={report['summary']['raw_labels']}; "
        f"proposed_primary_labels={report['summary']['proposed_primary_analysis_labels']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
