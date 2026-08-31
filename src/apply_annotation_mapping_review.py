from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import yaml

from .annotation_mapping import MAPPING_COLUMNS, file_sha256, load_label_mapping
from .utils import ensure_dir

FINAL_STATUSES = {"accept_as_proposed", "modify"}


def _load_review(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"Annotation review must be a mapping: {path}")
    return payload


def apply_review(
    *,
    review_path: Path,
    mapping_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    review = _load_review(review_path)
    metadata = review.get("review_metadata", {})
    if metadata.get("status") != "approved_for_application":
        raise RuntimeError("Annotation review is not approved for mapping application.")
    if not bool(metadata.get("outcome_blind_confirmation", False)):
        raise RuntimeError("Annotation review lacks outcome-blind confirmation.")
    if (
        not str(metadata.get("reviewer", "")).strip()
        or not str(metadata.get("review_date", "")).strip()
    ):
        raise RuntimeError("Annotation review requires reviewer identity and review date.")

    denominator = review.get("composition_denominator_decision", {})
    if denominator.get("status") != "approved" or denominator.get("selected_option") != (
        "all_qc_passed_cells_with_other_unresolved"
    ):
        raise RuntimeError("The approved composition-denominator decision is missing.")
    eligibility = review.get("method_eligibility_review", {})
    if eligibility.get("status") != "approved" or eligibility.get("decision") != (
        "accept_current_method_specific_contracts"
    ):
        raise RuntimeError("The approved method-eligibility decision is missing.")

    evidence = review.get("evidence_contract", {})
    expected_start_hash = str(evidence.get("mapping_sha256_at_review_start", ""))
    observed_start_hash = file_sha256(mapping_path)
    if not expected_start_hash or observed_start_hash != expected_start_hash:
        raise RuntimeError("Annotation mapping changed after the human review started.")

    mapping = load_label_mapping(mapping_path)
    decisions = review.get("decisions")
    if not isinstance(decisions, list) or not decisions:
        raise ValueError("Annotation review contains no label decisions.")
    raw_labels = [str(item.get("raw_label", "")).strip() for item in decisions]
    if len(raw_labels) != len(set(raw_labels)):
        raise ValueError("Annotation review contains duplicate raw-label decisions.")
    mapped_labels = set(mapping["raw_label"].astype(str))
    reviewed_labels = set(raw_labels)
    if mapped_labels != reviewed_labels:
        raise ValueError(
            "Annotation review does not match the mapping; "
            f"missing={sorted(mapped_labels - reviewed_labels)}; "
            f"extra={sorted(reviewed_labels - mapped_labels)}"
        )

    decision_lookup: dict[str, dict[str, str]] = {}
    for item in decisions:
        raw_label = str(item["raw_label"]).strip()
        final = item.get("final_decision", {})
        status = str(final.get("status", "")).strip()
        analysis_label = str(final.get("analysis_label", "")).strip()
        disposition = str(final.get("disposition", "")).strip()
        rationale = str(final.get("rationale", "")).strip()
        if status not in FINAL_STATUSES:
            raise RuntimeError(
                f"Annotation decision remains unresolved for {raw_label!r}: {status}"
            )
        if disposition not in {"primary", "exploratory", "unresolved"}:
            raise ValueError(f"Invalid final disposition for {raw_label!r}: {disposition}")
        if disposition != "unresolved" and not analysis_label:
            raise ValueError(f"Final analysis label is missing for {raw_label!r}.")
        if disposition == "unresolved" and analysis_label:
            raise ValueError(f"Unresolved label {raw_label!r} must not have an analysis label.")
        if not rationale:
            raise ValueError(f"Final rationale is missing for {raw_label!r}.")
        decision_lookup[raw_label] = {
            "analysis_label": analysis_label,
            "disposition": disposition,
            "rationale": rationale,
            "review_status": "approved",
        }

    applied = mapping.copy()
    for row_index, row in applied.iterrows():
        decision = decision_lookup[str(row["raw_label"])]
        for column, value in decision.items():
            applied.at[row_index, column] = value
    applied = applied[list(MAPPING_COLUMNS)]

    ensure_dir(output_path.parent)
    temporary = output_path.with_name(output_path.name + ".tmp")
    applied.to_csv(temporary, index=False, lineterminator="\n")
    validated = load_label_mapping(temporary, observed_labels=mapped_labels)
    if not validated["review_status"].eq("approved").all():
        raise RuntimeError("Applied annotation mapping is not fully approved.")
    os.replace(temporary, output_path)
    return {
        "rows": int(len(validated)),
        "mapping_version": str(validated["mapping_version"].iloc[0]),
        "mapping_sha256": file_sha256(output_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = apply_review(
        review_path=Path(args.review),
        mapping_path=Path(args.mapping),
        output_path=Path(args.out),
    )
    print(
        "[apply_annotation_mapping_review] applied "
        f"{result['rows']} approved rows; sha256={result['mapping_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
