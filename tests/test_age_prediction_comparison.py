from __future__ import annotations

import pandas as pd
import pytest

from src.age_prediction import _build_model_comparison_summary


def _predictions(visit_label: str, cell_type_label: str) -> pd.DataFrame:
    rows = []
    for model, visit_errors, cell_type_errors in [
        ("ridge", [4, 8, 2, 5], [12, 14, 16]),
        ("nested_selected", [1, 3, 4, 5 + 1e-10], [20, 20, 20]),
    ]:
        for level, folds, errors in [
            (visit_label, [0, 0, 1, 2], visit_errors),
            (cell_type_label, [0, 1, 2], cell_type_errors),
        ]:
            rows.extend(
                {
                    "model": model,
                    "params": "{}",
                    "evaluation_level": level,
                    "fold": fold,
                    "abs_error": error,
                }
                for fold, error in zip(folds, errors)
            )
    return pd.DataFrame(rows)


def _metrics(primary_level: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "model": ["ridge", "nested_selected"],
            "params": ["{}", "{}"],
            "primary_level": [primary_level, primary_level],
            "is_best": [False, True],
            "is_best_ridge": [True, False],
        }
    )


@pytest.mark.parametrize(
    ("visit_label", "cell_type_label", "primary_level", "expected_delta", "expected_counts"),
    [
        ("sample_unit", "sample_unit_celltype", "sample_unit", 2 / 3, (1, 1, 1)),
        ("donor", "donor_celltype", "donor", 2 / 3, (1, 1, 1)),
        ("sample_unit", "sample_unit_celltype", "donor", 2 / 3, (1, 1, 1)),
        ("sample_unit", "sample_unit_celltype", "sample_unit_celltype", -6, (0, 3, 0)),
        ("donor", "donor_celltype", "donor_celltype", -6, (0, 3, 0)),
    ],
)
def test_model_comparison_uses_one_level_and_one_outcome_per_fold(
    visit_label: str,
    cell_type_label: str,
    primary_level: str,
    expected_delta: float,
    expected_counts: tuple[int, int, int],
) -> None:
    summary = _build_model_comparison_summary(
        _predictions(visit_label, cell_type_label), _metrics(primary_level), {}
    ).set_index("candidate_model")

    nested = summary.loc["nested_selected"]
    assert nested["delta_mae_mean"] == pytest.approx(expected_delta)
    assert nested["n_folds_compared"] == 3
    assert (
        nested["winner_count_candidate"],
        nested["winner_count_reference"],
        nested["tie_count"],
    ) == expected_counts
    ridge = summary.loc["ridge"]
    assert ridge["n_folds_compared"] == ridge["tie_count"] == 3
    assert ridge["delta_mae_boot_ci_low"] == ridge["delta_mae_boot_ci_high"] == 0


def test_model_comparison_refuses_other_levels_when_primary_predictions_are_missing() -> None:
    predictions = _predictions("sample_unit", "sample_unit_celltype")
    cell_types = predictions.loc[predictions["evaluation_level"].eq("sample_unit_celltype")]

    with pytest.raises(ValueError, match="exactly one set of sample_unit predictions"):
        _build_model_comparison_summary(cell_types, _metrics("sample_unit"), {})


def test_model_comparison_refuses_incomplete_pairs() -> None:
    predictions = _predictions("sample_unit", "sample_unit_celltype")
    missing = (
        predictions["model"].eq("nested_selected")
        & predictions["fold"].eq(2)
        & predictions["evaluation_level"].eq("sample_unit")
    )

    with pytest.raises(ValueError, match="same complete folds"):
        _build_model_comparison_summary(predictions.loc[~missing], _metrics("sample_unit"), {})
