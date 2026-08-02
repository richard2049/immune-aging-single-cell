from __future__ import annotations

import unittest

from src.draft_corrected_public_assets import (
    _validate_composition_selection,
)


class DraftCorrectedPublicAssetTests(unittest.TestCase):
    def test_composition_selection_requires_retained_signals(self) -> None:
        signals = [
            "composition::A",
            "composition::B",
            "composition::C",
        ]
        decisions = {signal: {"disposition": "retain"} for signal in signals}

        observed = _validate_composition_selection(signals, decisions)

        self.assertEqual(observed, ["A", "B", "C"])

    def test_composition_selection_rejects_excluded_signal(self) -> None:
        signals = [
            "composition::A",
            "composition::B",
            "composition::C",
        ]
        decisions = {
            "composition::A": {"disposition": "retain"},
            "composition::B": {"disposition": "exclude"},
            "composition::C": {"disposition": "retain"},
        }

        with self.assertRaisesRegex(ValueError, "non-retained"):
            _validate_composition_selection(signals, decisions)


if __name__ == "__main__":
    unittest.main()
