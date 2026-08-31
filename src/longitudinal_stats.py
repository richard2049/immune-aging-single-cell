from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .scientific_guardrails import NonEstimableDesignError, build_complete_nuisance_design


class GEENonConvergenceError(RuntimeError):
    """Raised when a GEE fit exhausts its prespecified iteration budget."""


def select_one_sample_per_subject(
    frame: pd.DataFrame,
    *,
    subject_col: str,
    sample_col: str,
    age_col: str = "age",
    rule: str = "earliest_age_then_sample_id",
) -> set[str]:
    """Return a deterministic sample-unit set with at most one sample per subject."""
    if rule != "earliest_age_then_sample_id":
        raise ValueError(f"Unsupported one-sample-per-subject rule: {rule!r}.")
    required = [subject_col, sample_col, age_col]
    missing = [column for column in required if column not in frame]
    if missing:
        raise KeyError(f"One-sample sensitivity is missing columns: {missing}.")

    metadata = frame[required].copy()
    if metadata.isna().any(axis=None):
        raise ValueError("One-sample sensitivity metadata contain missing values.")
    metadata[subject_col] = metadata[subject_col].astype("string").str.strip()
    metadata[sample_col] = metadata[sample_col].astype("string").str.strip()
    metadata[age_col] = pd.to_numeric(metadata[age_col], errors="coerce")
    if (
        metadata[subject_col].eq("").any()
        or metadata[sample_col].eq("").any()
        or not np.isfinite(metadata[age_col].to_numpy(dtype=float)).all()
    ):
        raise ValueError("One-sample sensitivity metadata contain blank or invalid values.")

    sample_contract = metadata.groupby(sample_col, observed=False).agg(
        n_subjects=(subject_col, "nunique"),
        n_ages=(age_col, "nunique"),
    )
    if bool(sample_contract["n_subjects"].ne(1).any()) or bool(
        sample_contract["n_ages"].ne(1).any()
    ):
        raise ValueError("A sample unit maps to conflicting subject or age values.")

    unique_samples = metadata.drop_duplicates(sample_col).sort_values(
        [subject_col, age_col, sample_col],
        kind="stable",
    )
    selected = unique_samples.drop_duplicates(subject_col, keep="first")
    return set(selected[sample_col].astype(str))


def fit_independent_age_model(
    frame: pd.DataFrame,
    *,
    outcome_col: str,
    covariate_cols: list[str],
    family: str,
    weights_col: str | None = None,
    context: str,
) -> dict[str, Any]:
    """Fit the pre-specified one-sample-per-subject sensitivity model."""
    required = [outcome_col, "age", *covariate_cols]
    if weights_col:
        required.append(weights_col)
    missing_columns = [column for column in required if column not in frame]
    if missing_columns:
        raise KeyError(f"{context}: missing model columns: {missing_columns}.")
    data = frame[required].copy()
    if data.isna().any(axis=None):
        raise ValueError(f"{context}: model inputs contain missing values.")

    age = pd.to_numeric(data["age"], errors="coerce").to_numpy(dtype=float)
    outcome = pd.to_numeric(data[outcome_col], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(age).all() or np.ptp(age) <= 0:
        raise ValueError(f"{context}: age must be finite and variable.")
    if not np.isfinite(outcome).all():
        raise ValueError(f"{context}: outcome contains non-finite values.")

    nuisance = (
        build_complete_nuisance_design(data[covariate_cols], context=f"{context} fixed effects")
        if covariate_cols
        else np.ones((len(data), 1), dtype=float)
    )
    age_decade = ((age - float(np.mean(age))) / 10.0).reshape(-1, 1)
    design = np.column_stack([nuisance, age_decade])
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise NonEstimableDesignError(
            f"{context}: fixed-effect design is rank deficient "
            f"(rank={rank}, columns={design.shape[1]})."
        )

    weights = None
    if weights_col:
        weights = pd.to_numeric(data[weights_col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(weights).all() or np.any(weights <= 0):
            raise ValueError(f"{context}: model weights must be finite and positive.")

    import statsmodels.api as sm

    family_key = str(family).strip().lower()
    if family_key == "binomial":
        if np.any((outcome < 0) | (outcome > 1)):
            raise ValueError(f"{context}: binomial outcome must lie in [0, 1].")
        model_family = sm.families.Binomial()
    elif family_key == "gaussian":
        model_family = sm.families.Gaussian()
    else:
        raise ValueError(f"{context}: unsupported model family {family!r}.")

    model = sm.GLM(
        endog=outcome,
        exog=design,
        family=model_family,
        freq_weights=weights,
    )
    fit = model.fit(maxiter=200, cov_type="HC0")
    if not bool(getattr(fit, "converged", False)):
        raise RuntimeError(f"{context}: sensitivity model did not converge.")

    age_index = design.shape[1] - 1
    ci = np.asarray(fit.conf_int(), dtype=float)[age_index]
    return {
        "model": f"GLM-{family_key}-HC0",
        "effect_per_10y": float(fit.params[age_index]),
        "effect_per_10y_ci_low": float(ci[0]),
        "effect_per_10y_ci_high": float(ci[1]),
        "pvalue": float(fit.pvalues[age_index]),
        "n_observations": int(len(data)),
        "design_columns": int(design.shape[1]),
        "design_rank": rank,
        "converged": True,
    }


def fit_subject_aware_gee(
    frame: pd.DataFrame,
    *,
    outcome_col: str,
    subject_col: str,
    covariate_cols: list[str],
    family: str,
    weights_col: str | None = None,
    working_correlation: str = "exchangeable",
    context: str,
) -> dict[str, Any]:
    """Fit a population-average age association with subject-clustered GEE."""
    required = [outcome_col, subject_col, "age", *covariate_cols]
    if weights_col:
        required.append(weights_col)
    missing_columns = [column for column in required if column not in frame]
    if missing_columns:
        raise KeyError(f"{context}: missing GEE columns: {missing_columns}.")

    data = frame[required].copy()
    missing_values = data.isna().sum()
    missing_values = missing_values[missing_values.gt(0)]
    if not missing_values.empty:
        detail = ", ".join(f"{column}={int(count)}" for column, count in missing_values.items())
        raise ValueError(f"{context}: GEE inputs contain missing values ({detail}).")

    subject = data[subject_col].astype("string").str.strip()
    if subject.eq("").any():
        raise ValueError(f"{context}: subject identifiers contain blank values.")
    n_subjects = int(subject.nunique())
    if n_subjects < 3:
        raise ValueError(f"{context}: at least three subjects are required for GEE.")

    age = pd.to_numeric(data["age"], errors="coerce").to_numpy(dtype=float)
    outcome = pd.to_numeric(data[outcome_col], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(age).all() or np.ptp(age) <= 0:
        raise ValueError(f"{context}: age must be finite and variable.")
    if not np.isfinite(outcome).all():
        raise ValueError(f"{context}: outcome contains non-finite values.")

    if covariate_cols:
        nuisance = build_complete_nuisance_design(
            data[covariate_cols], context=f"{context} fixed effects"
        )
    else:
        nuisance = np.ones((len(data), 1), dtype=float)
    age_decade = ((age - float(np.mean(age))) / 10.0).reshape(-1, 1)
    design = np.column_stack([nuisance, age_decade])
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise NonEstimableDesignError(
            f"{context}: fixed-effect design is rank deficient "
            f"(rank={rank}, columns={design.shape[1]})."
        )

    weights = None
    if weights_col:
        weights = pd.to_numeric(data[weights_col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(weights).all() or np.any(weights <= 0):
            raise ValueError(f"{context}: GEE weights must be finite and positive.")

    import statsmodels.api as sm

    family_key = str(family).strip().lower()
    if family_key == "binomial":
        if np.any((outcome < 0) | (outcome > 1)):
            raise ValueError(f"{context}: binomial outcome must lie in [0, 1].")
        model_family = sm.families.Binomial()
    elif family_key == "gaussian":
        model_family = sm.families.Gaussian()
    else:
        raise ValueError(f"{context}: unsupported GEE family {family!r}.")

    correlation_key = str(working_correlation).strip().lower()
    if correlation_key == "exchangeable":
        covariance_structure = sm.cov_struct.Exchangeable()
    elif correlation_key == "independence":
        covariance_structure = sm.cov_struct.Independence()
    else:
        raise ValueError(f"{context}: unsupported GEE working correlation {working_correlation!r}.")

    model = sm.GEE(
        endog=outcome,
        exog=design,
        groups=subject.to_numpy(dtype=str),
        family=model_family,
        cov_struct=covariance_structure,
        weights=weights,
    )
    fit = model.fit(maxiter=200, cov_type="robust")
    if not bool(getattr(fit, "converged", False)):
        raise GEENonConvergenceError(
            f"{context}: GEE with {correlation_key} working correlation did not converge."
        )

    age_index = design.shape[1] - 1
    ci = np.asarray(fit.conf_int(), dtype=float)[age_index]
    inferential_values = np.asarray(
        [fit.params[age_index], ci[0], ci[1], fit.pvalues[age_index], fit.bse[age_index]],
        dtype=float,
    )
    if not np.isfinite(inferential_values).all():
        raise RuntimeError(
            f"{context}: GEE with {correlation_key} working correlation produced "
            "non-finite inference."
        )
    dependence_parameters = getattr(fit.cov_struct, "dep_params", None)
    dependence = (
        np.asarray(dependence_parameters, dtype=float).reshape(-1)
        if dependence_parameters is not None
        else np.asarray([], dtype=float)
    )
    if correlation_key == "exchangeable" and (
        dependence.size != 1 or not np.isfinite(dependence[0])
    ):
        raise RuntimeError(f"{context}: exchangeable GEE produced an invalid working correlation.")
    return {
        "model": f"GEE-{family_key}-{correlation_key}",
        "effect_per_10y": float(fit.params[age_index]),
        "effect_per_10y_ci_low": float(ci[0]),
        "effect_per_10y_ci_high": float(ci[1]),
        "pvalue": float(fit.pvalues[age_index]),
        "standard_error": float(fit.bse[age_index]),
        "n_observations": int(len(data)),
        "n_subjects": n_subjects,
        "design_columns": int(design.shape[1]),
        "design_rank": rank,
        "working_correlation_structure": correlation_key,
        "working_correlation": (float(dependence[0]) if dependence.size else np.nan),
        "covariance_type": "robust_subject_clustered",
        "converged": True,
        "effect_scale": (
            "log_odds_per_10_years" if family_key == "binomial" else "outcome_units_per_10_years"
        ),
    }
