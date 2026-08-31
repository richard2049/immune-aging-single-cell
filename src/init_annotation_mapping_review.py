from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .annotation_mapping import file_sha256, load_label_mapping
from .utils import ensure_dir


def initialize_review(
    *,
    mapping_path: Path,
    review_table_path: Path,
    qualification_report_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    if output_path.exists():
        raise FileExistsError(
            f"Review record already exists and will not be overwritten: {output_path}"
        )
    mapping = load_label_mapping(mapping_path)
    review_table = pd.read_csv(review_table_path, dtype={"raw_label": "string"})
    if review_table["raw_label"].duplicated().any():
        raise ValueError("Annotation review table contains duplicate raw labels.")
    mapped = set(mapping["raw_label"].astype(str))
    observed = set(review_table["raw_label"].dropna().astype(str))
    if mapped != observed:
        raise ValueError(
            "Annotation review table does not match the mapping; "
            f"missing={sorted(mapped - observed)}; extra={sorted(observed - mapped)}"
        )

    disposition_order = {"primary": 0, "exploratory": 1, "unresolved": 2}
    ordered = mapping.assign(
        _order=mapping["disposition"].map(disposition_order).astype(int)
    ).sort_values(["_order", "analysis_label", "raw_label"])
    decisions: list[dict[str, Any]] = []
    for review_order, row in enumerate(ordered.itertuples(index=False), start=1):
        decisions.append(
            {
                "review_order": review_order,
                "raw_label": str(row.raw_label),
                "proposed_analysis_label": str(row.analysis_label),
                "proposed_disposition": str(row.disposition),
                "mapping_rationale": str(row.rationale),
                "evidence_row": str(row.raw_label),
                "human_assessment": {
                    "marker_coherence": "not_started",
                    "cluster_coherence": "not_started",
                    "sample_subject_support": "not_started",
                    "technical_distribution": "not_started",
                    "pbmc_plausibility": "not_started",
                },
                "final_decision": {
                    "status": "pending_review",
                    "analysis_label": str(row.analysis_label),
                    "disposition": str(row.disposition),
                    "rationale": "",
                },
            }
        )

    payload: dict[str, Any] = {
        "review_metadata": {
            "status": "in_progress",
            "reviewer": "",
            "review_date": "",
            "outcome_blind_confirmation": False,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "mapping_version": str(mapping["mapping_version"].iloc[0]),
        },
        "evidence_contract": {
            "mapping_path": mapping_path.as_posix(),
            "mapping_sha256_at_review_start": file_sha256(mapping_path),
            "review_table_path": review_table_path.as_posix(),
            "review_table_sha256": file_sha256(review_table_path),
            "qualification_report_path": qualification_report_path.as_posix(),
            "qualification_report_sha256": file_sha256(qualification_report_path),
            "forbidden_evidence": [
                "age-association estimates",
                "effect directions",
                "p-values or FDR",
                "prediction performance",
                "gene-level differential-expression results",
            ],
        },
        "composition_denominator_decision": {
            "status": "pending_review",
            "selected_option": "",
            "recommended_option": "all_qc_passed_cells_with_other_unresolved",
            "allowed_options": [
                "all_qc_passed_cells_with_other_unresolved",
                "primary_mapped_cells_only_with_explicit_label",
            ],
            "rationale": "",
        },
        "method_eligibility_review": {
            "status": "pending_review",
            "decision": "",
            "recommended_decision": "accept_current_method_specific_contracts",
            "rationale": "",
        },
        "decisions": decisions,
        "final_approval": {
            "status": "pending_review",
            "mapping_applied": False,
            "mapping_sha256_after_application": "",
            "reviewer_confirmation": "",
        },
    }
    ensure_dir(output_path.parent)
    output_path.write_text(
        yaml.safe_dump(payload, sort_keys=False, allow_unicode=False, width=100),
        encoding="utf-8",
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--review-table", required=True)
    parser.add_argument("--qualification-report", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    payload = initialize_review(
        mapping_path=Path(args.mapping),
        review_table_path=Path(args.review_table),
        qualification_report_path=Path(args.qualification_report),
        output_path=Path(args.out),
    )
    print(
        "[init_annotation_mapping_review] created "
        f"{len(payload['decisions'])} pending review entries",
        flush=True,
    )


if __name__ == "__main__":
    main()
