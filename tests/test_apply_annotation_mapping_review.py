from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from src.annotation_mapping import file_sha256, load_label_mapping
from src.apply_annotation_mapping_review import apply_review


class ApplyAnnotationMappingReviewTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> tempfile.TemporaryDirectory:
        root = Path(".tmp/tests")
        root.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=root)

    @staticmethod
    def _write_fixture(root: Path) -> tuple[Path, Path]:
        mapping_path = root / "mapping.csv"
        pd.DataFrame(
            [
                {
                    "raw_label": "A",
                    "analysis_label": "Old T label",
                    "disposition": "primary",
                    "rationale": "Initial proposal A",
                    "review_status": "pending_review",
                    "mapping_version": "fixture-v1",
                },
                {
                    "raw_label": "B",
                    "analysis_label": "Old rare label",
                    "disposition": "exploratory",
                    "rationale": "Initial proposal B",
                    "review_status": "pending_review",
                    "mapping_version": "fixture-v1",
                },
            ]
        ).to_csv(mapping_path, index=False)
        review_path = root / "review.yml"
        review = {
            "review_metadata": {
                "status": "approved_for_application",
                "reviewer": "reviewer",
                "review_date": "2026-08-25",
                "outcome_blind_confirmation": True,
            },
            "evidence_contract": {
                "mapping_sha256_at_review_start": file_sha256(mapping_path),
            },
            "composition_denominator_decision": {
                "status": "approved",
                "selected_option": "all_qc_passed_cells_with_other_unresolved",
            },
            "method_eligibility_review": {
                "status": "approved",
                "decision": "accept_current_method_specific_contracts",
            },
            "decisions": [
                {
                    "raw_label": "A",
                    "final_decision": {
                        "status": "modify",
                        "analysis_label": "T cells",
                        "disposition": "primary",
                        "rationale": "Reviewed broad T lineage",
                    },
                },
                {
                    "raw_label": "B",
                    "final_decision": {
                        "status": "accept_as_proposed",
                        "analysis_label": "",
                        "disposition": "unresolved",
                        "rationale": "Reviewed exclusion from primary labels",
                    },
                },
            ],
        }
        review_path.write_text(yaml.safe_dump(review, sort_keys=False), encoding="utf-8")
        return mapping_path, review_path

    def test_review_application_preserves_raw_labels_and_approves_exclusion(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            mapping_path, review_path = self._write_fixture(root)
            output_path = root / "approved.csv"

            result = apply_review(
                review_path=review_path,
                mapping_path=mapping_path,
                output_path=output_path,
            )

            approved = load_label_mapping(output_path)
            self.assertEqual(approved["raw_label"].tolist(), ["A", "B"])
            self.assertTrue(approved["review_status"].eq("approved").all())
            unresolved = approved.loc[approved["raw_label"].eq("B")].iloc[0]
            self.assertEqual(unresolved["analysis_label"], "")
            self.assertEqual(unresolved["disposition"], "unresolved")
            self.assertEqual(result["mapping_sha256"], file_sha256(output_path))

    def test_review_application_rejects_mapping_changed_after_review(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            mapping_path, review_path = self._write_fixture(root)
            mapping = pd.read_csv(mapping_path)
            mapping.loc[0, "rationale"] = "Changed after review"
            mapping.to_csv(mapping_path, index=False)

            with self.assertRaisesRegex(RuntimeError, "changed after"):
                apply_review(
                    review_path=review_path,
                    mapping_path=mapping_path,
                    output_path=root / "approved.csv",
                )


if __name__ == "__main__":
    unittest.main()
