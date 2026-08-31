from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from src.scientific_guardrails import NonEstimableDesignError
from src.sensitivity_age import _outputs_complete
from src.signature_age import _association_table
from src.supplementary_age_plots import _prepare_comp_effects


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
