from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import pandas as pd

from src.review_corrected_results import (
    _load_comparison,
    apply_human_dispositions,
    build_promotion_manifest,
    classify_association,
)


CRITERIA = {
    "require_effect_ci_excludes_zero": True,
    "require_all_sensitivity_scenarios": True,
    "require_direction_consistency": True,
    "require_fdr_support_all_scenarios": True,
}


def _row(**overrides: bool) -> dict[str, bool]:
    row = {
        "corrected_fdr_pass": True,
        "min_support_pass": True,
        "effect_ci_excludes_zero": True,
        "grouping_valid": True,
        "sensitivity_complete_pass": True,
        "direction_consistency_pass": True,
        "all_sensitivity_fdr_pass": True,
        "sensitivity_grouping_valid": True,
    }
    row.update(overrides)
    return row


class CorrectedResultReviewTests(unittest.TestCase):
    def test_complete_stable_evidence_is_candidate_for_human_review(
        self,
    ) -> None:
        observed = classify_association(_row(), CRITERIA)

        self.assertEqual(observed, "candidate_for_human_review")

    def test_sensitivity_instability_remains_exploratory(self) -> None:
        observed = classify_association(
            _row(direction_consistency_pass=False),
            CRITERIA,
        )

        self.assertEqual(observed, "exploratory_sensitivity_limited")

    def test_failed_corrected_fdr_is_not_prioritized(self) -> None:
        observed = classify_association(
            _row(corrected_fdr_pass=False),
            CRITERIA,
        )

        self.assertEqual(observed, "not_prioritized")

    def test_incorrect_grouping_is_not_prioritized(self) -> None:
        observed = classify_association(
            _row(grouping_valid=False),
            CRITERIA,
        )

        self.assertEqual(observed, "not_prioritized")

    def test_comparison_normalizes_empty_signature(self) -> None:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=root) as tmp:
            path = Path(tmp) / "comparison.csv"
            pd.DataFrame(
                {
                    "analysis": ["composition"],
                    "cell_type": ["Test cells"],
                    "signature": [None],
                    "comparison_status": ["stable_direction"],
                }
            ).to_csv(path, index=False)

            observed = _load_comparison(path)

        self.assertIn(("composition", "Test cells", ""), observed)

    def test_human_disposition_is_applied_without_approving_claim(
        self,
    ) -> None:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        evidence = pd.DataFrame(
            {
                "signal_id": ["composition::Test cells"],
                "automated_screen_status": [
                    "candidate_for_human_review"
                ],
                "human_disposition": ["pending"],
                "public_claim_approved": [False],
            }
        )
        with tempfile.TemporaryDirectory(dir=root) as tmp:
            path = Path(tmp) / "dispositions.yml"
            path.write_text(
                """
review_metadata:
  status: dispositions_approved_public_wording_pending
  reviewer: Repository owner
  review_date: "2026-07-30"
  public_scope: conservative_core
  public_claims_approved: false
  asset_replacements_approved: false
decisions:
  - signal_id: "composition::Test cells"
    disposition: retain
    presentation_tier: core
    rationale: Stable evidence.
""".lstrip(),
                encoding="utf-8",
            )

            observed, metadata = apply_human_dispositions(evidence, path)

        self.assertEqual(observed.loc[0, "human_disposition"], "retain")
        self.assertEqual(observed.loc[0, "presentation_tier"], "core")
        self.assertFalse(observed.loc[0, "public_claim_approved"])
        self.assertTrue(metadata["candidate_dispositions_complete"])

    def test_human_disposition_rejects_unknown_signal(self) -> None:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        evidence = pd.DataFrame(
            {
                "signal_id": ["composition::Known cells"],
                "automated_screen_status": [
                    "candidate_for_human_review"
                ],
                "human_disposition": ["pending"],
            }
        )
        with tempfile.TemporaryDirectory(dir=root) as tmp:
            path = Path(tmp) / "dispositions.yml"
            path.write_text(
                """
review_metadata:
  status: dispositions_approved_public_wording_pending
  reviewer: Repository owner
  review_date: "2026-07-30"
  public_scope: conservative_core
  public_claims_approved: false
  asset_replacements_approved: false
decisions:
  - signal_id: "composition::Unknown cells"
    disposition: retain
    presentation_tier: core
    rationale: Stable evidence.
""".lstrip(),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "must match one evidence row",
            ):
                apply_human_dispositions(evidence, path)

    def test_promotion_manifest_keeps_public_changes_blocked(
        self,
    ) -> None:
        evidence = pd.DataFrame(
            {
                "signal_id": ["composition::Test cells"],
                "automated_screen_status": [
                    "candidate_for_human_review"
                ],
                "human_disposition": ["retain"],
            }
        )
        summary = {
            "review_status": "draft_asset_design_approved",
            "review_metadata": {
                "reviewer": "Repository owner",
                "review_date": "2026-07-30",
                "public_scope": "conservative_core",
                "draft_asset_decisions": {
                    "signature_figure": "omit_from_curated_assets"
                },
            },
        }

        observed = build_promotion_manifest(
            evidence,
            Path("results/corrected"),
            summary,
        )

        self.assertFalse(observed["public_changes_authorized"])
        signature = next(
            item
            for item in observed["assets"]
            if item["role"] == "signature"
        )
        self.assertEqual(
            signature["status"],
            "not_selected_for_curated_assets",
        )


if __name__ == "__main__":
    unittest.main()
