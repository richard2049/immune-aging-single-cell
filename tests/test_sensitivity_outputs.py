from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from scipy.special import logit

from src.composition_age import _compute_trends, _fit_one_sample_composition
from src.scientific_guardrails import NonEstimableDesignError
from src.sensitivity_age import _outputs_complete
from src.signature_age import _association_table
from src.supplementary_age_plots import _composition_effect_xlabel, _prepare_comp_effects


def test_non_estimable_design_error_remains_a_value_error() -> None:
    assert issubclass(NonEstimableDesignError, ValueError)


def _signature_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_type": ["B cells"] * 4,
            "donor_id": ["T1", "T2", "T3", "T4"],
            "subject_id": ["S1", "S2", "S3", "S4"],
            "age": [20.0, 30.0, 40.0, 50.0],
            "score__ifn_response": [0.1, 0.3, 0.2, 0.6],
        }
    )


def _primary_gee_result() -> dict[str, object]:
    return {
        "pvalue": 0.2,
        "effect_per_10y": 0.1,
        "effect_per_10y_ci_low": -0.1,
        "effect_per_10y_ci_high": 0.3,
        "model": "GEE-gaussian-exchangeable",
        "effect_scale": "score_units_per_10_years",
        "standard_error": 0.1,
        "working_correlation_structure": "exchangeable",
        "working_correlation": 0.2,
        "covariance_type": "robust_subject_clustered",
        "fallback_used": False,
        "fallback_reason": "",
        "fallback_requires_manual_review": False,
    }


def test_signature_secondary_non_estimability_is_recorded() -> None:
    with (
        patch("src.signature_age._fit_signature_gee", return_value=_primary_gee_result()),
        patch(
            "src.signature_age.fit_independent_age_model",
            side_effect=NonEstimableDesignError("fixture: rank deficient"),
        ),
    ):
        observed = _association_table(
            _signature_fixture(),
            score_cols=["score__ifn_response"],
            min_donors_per_celltype=3,
            min_age_span=20.0,
            covariate_cols=[],
            adjust_covariates=False,
            bootstrap_iterations=10,
            bootstrap_ci=0.95,
            seed=42,
            subject_aware=True,
        )

    assert observed.loc[0, "sensitivity_status"] == "non_estimable_rank_deficient"
    assert "rank deficient" in observed.loc[0, "sensitivity_reason"]
    assert pd.isna(observed.loc[0, "sensitivity_effect_per_10y"])
    assert observed.loc[0, "association_model"] == "GEE-gaussian-exchangeable"


def test_signature_secondary_does_not_mask_unrelated_errors() -> None:
    with (
        patch("src.signature_age._fit_signature_gee", return_value=_primary_gee_result()),
        patch(
            "src.signature_age.fit_independent_age_model",
            side_effect=RuntimeError("unexpected fit failure"),
        ),
        pytest.raises(RuntimeError, match="unexpected fit failure"),
    ):
        _association_table(
            _signature_fixture(),
            score_cols=["score__ifn_response"],
            min_donors_per_celltype=3,
            min_age_span=20.0,
            covariate_cols=[],
            adjust_covariates=False,
            bootstrap_iterations=10,
            bootstrap_ci=0.95,
            seed=42,
            subject_aware=True,
        )


def test_outputs_complete_requires_every_nonempty_file(tmp_path: Path) -> None:
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    first.write_text("value\n1\n", encoding="utf-8")
    assert not _outputs_complete([first, second])
    second.touch()
    assert not _outputs_complete([first, second])
    second.write_text("value\n2\n", encoding="utf-8")
    assert _outputs_complete([first, second])


@pytest.fixture
def composition_people() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "subject_id": [f"S{index}" for index in range(8)],
            "age": [30.0] * 4 + [70.0] * 4,
            "fraction": [0.1, 0.15, 0.2, 0.25, 0.12, 0.2, 0.22, 0.3],
            "total_cells": [100, 200, 300, 400, 200, 100, 400, 300],
        }
    )


def test_composition_sensitivity_uses_person_level_uncertainty(
    composition_people: pd.DataFrame,
) -> None:
    people = composition_people
    # With two age groups the model is saturated. Each group logit is determined
    # by its weighted mean, and its robust variance by the people in the group.
    group_logits, group_variances = [], []
    for _, group in people.groupby("age"):
        weights = group["total_cells"].to_numpy()
        fractions = group["fraction"].to_numpy()
        mean = np.average(fractions, weights=weights)
        group_logits.append(logit(mean))
        group_variances.append(
            np.sum((weights * (fractions - mean)) ** 2) / (weights.sum() * mean * (1 - mean)) ** 2
        )
    expected_effect = (group_logits[1] - group_logits[0]) / 4
    expected_se = np.sqrt(sum(group_variances)) / 4
    frame = people.assign(
        donor_id=people["subject_id"],
        cell_type="B cells",
        n_cells=(people["fraction"] * people["total_cells"]).round().astype(int),
    )
    row = _compute_trends(
        frame,
        min_donors_per_celltype=3,
        covariate_cols=[],
        adjust_covariates=False,
        bootstrap_iterations=10,
        bootstrap_ci=0.95,
        seed=42,
        subject_aware=True,
    ).iloc[0]
    assert row["association_model"] == "GEE-binomial-exchangeable"
    assert row["sensitivity_model"] == "GEE-binomial-independence"
    assert row["sensitivity_covariance_type"] == "robust_subject_clustered"
    assert row["sensitivity_n_subjects"] == len(people)
    assert row["sensitivity_effect_per_10y"] == pytest.approx(expected_effect, abs=1e-10)
    assert row["sensitivity_standard_error"] == pytest.approx(expected_se, abs=1e-10)

    scaled = _fit_one_sample_composition(
        people.assign(total_cells=people["total_cells"] * 100),
        covariate_cols=[],
        context="test scaled weights",
    )
    assert scaled["standard_error"] == pytest.approx(expected_se, abs=1e-10)


@pytest.mark.parametrize("subject_id", ["S1", "", None])
def test_composition_sensitivity_rejects_invalid_subject_groups(
    composition_people: pd.DataFrame, subject_id: str | None
) -> None:
    people = composition_people.copy()
    people.loc[0, "subject_id"] = subject_id
    with pytest.raises(ValueError, match="one observation per subject"):
        _fit_one_sample_composition(people, covariate_cols=[], context="test subject groups")


def test_composition_supplement_uses_primary_association_fields() -> None:
    source = pd.DataFrame(
        {
            "cell_type": ["B cells"],
            "slope_per_10y": [99.0],
            "slope_per_10y_ci_low": [98.0],
            "slope_per_10y_ci_high": [100.0],
            "spearman_fdr": [0.001],
            "association_effect_per_10y": [0.02],
            "association_effect_ci_low": [0.01],
            "association_effect_ci_high": [0.03],
            "association_fdr": [0.2],
            "fdr_significant": [False],
        }
    )

    observed = _prepare_comp_effects(source).iloc[0]

    assert observed["effect"] == 0.02
    assert observed["ci_low"] == 0.01
    assert observed["ci_high"] == 0.03
    assert observed["fdr"] == 0.2
    assert not bool(observed["significant"])


@pytest.mark.parametrize(
    ("scale", "expected"),
    [("log_odds_per_10_years", "log-odds"), ("fraction_per_10_years", "fraction units")],
)
def test_composition_forest_label_matches_recorded_scale(scale: str, expected: str) -> None:
    frame = pd.DataFrame({"effect_scale": [scale]})
    assert expected in _composition_effect_xlabel(frame)


@pytest.mark.parametrize(
    "scales", [["unspecified"], ["log_odds_per_10_years", "fraction_per_10_years"]]
)
def test_composition_forest_rejects_unknown_or_mixed_units(scales: list[str]) -> None:
    with pytest.raises(ValueError, match="single supported effect scale"):
        _composition_effect_xlabel(pd.DataFrame({"effect_scale": scales}))
