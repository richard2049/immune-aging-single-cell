from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype


def placeholders_allowed(config: dict) -> bool:
    """Return whether the active profile explicitly permits placeholder outputs."""
    return bool(config.get("run", {}).get("allow_placeholder_outputs", False))


def require_placeholder_permission(config: dict, reason: str) -> None:
    """Fail real-data profiles instead of converting analytical errors to files."""
    if not placeholders_allowed(config):
        raise RuntimeError(f"{reason} Placeholder outputs are disabled for this profile.")


def resolve_covariates(
    columns: Sequence[str],
    section: dict,
    *,
    section_name: str,
    defaults: Sequence[str] = (),
) -> list[str]:
    """Resolve covariates while enforcing explicitly configured columns."""
    configured = section.get("covariate_cols")
    explicit = configured is not None
    if explicit:
        if not isinstance(configured, list):
            raise TypeError(f"{section_name}.covariate_cols must be a list.")
        requested = [str(value) for value in configured]
    else:
        requested = [str(value) for value in defaults]

    duplicates = sorted({value for value in requested if requested.count(value) > 1})
    if duplicates:
        raise ValueError(f"{section_name}.covariate_cols contains duplicates: {duplicates}")

    available = set(columns)
    missing = [value for value in requested if value not in available]
    if explicit and missing:
        raise KeyError(
            f"Configured {section_name} covariates are missing: {missing}. "
            "Fix the metadata or approve a different model explicitly."
        )
    return [value for value in requested if value in available]


def validate_replicate_covariates(
    frame: pd.DataFrame,
    *,
    replicate_col: str,
    covariate_cols: Sequence[str],
    context: str,
) -> None:
    """Require complete, internally consistent replicate-level covariates."""
    if not covariate_cols:
        return
    missing = frame[list(covariate_cols)].isna().sum()
    missing = missing[missing.gt(0)]
    if not missing.empty:
        detail = ", ".join(f"{name}={int(count)}" for name, count in missing.items())
        raise ValueError(f"{context}: configured covariates contain missing values ({detail}).")

    for covariate in covariate_cols:
        levels = frame.groupby(replicate_col, observed=False)[covariate].nunique(dropna=False)
        inconsistent = levels[levels.gt(1)]
        if not inconsistent.empty:
            examples = inconsistent.index.astype(str).tolist()[:5]
            raise ValueError(
                f"{context}: covariate {covariate!r} is inconsistent within "
                f"biological replicates; examples: {examples}."
            )


def build_complete_nuisance_design(
    frame: pd.DataFrame,
    *,
    context: str,
) -> np.ndarray:
    """Build a complete, full-rank nuisance design with an intercept."""
    if frame.empty:
        raise ValueError(f"{context}: nuisance-covariate table is empty.")
    missing = frame.isna().sum()
    missing = missing[missing.gt(0)]
    if not missing.empty:
        detail = ", ".join(f"{name}={int(count)}" for name, count in missing.items())
        raise ValueError(
            f"{context}: configured covariates contain missing values "
            f"({detail}); silent imputation is not permitted."
        )

    matrices: list[np.ndarray] = []
    for column in frame.columns:
        values = frame[column]
        if values.nunique(dropna=False) < 2:
            raise ValueError(f"{context}: configured covariate {column!r} has no variation.")
        numeric = pd.to_numeric(values, errors="coerce")
        if is_numeric_dtype(values.dtype):
            array = numeric.to_numpy(dtype=float).reshape(-1, 1)
            if not np.isfinite(array).all():
                raise ValueError(f"{context}: numeric covariate {column!r} is not finite.")
            matrices.append(array)
            continue

        categories = values.astype("string")
        dummies = pd.get_dummies(
            categories,
            prefix=str(column),
            drop_first=True,
            dtype=float,
        )
        if dummies.shape[1]:
            matrices.append(dummies.to_numpy(dtype=float))

    intercept = np.ones((len(frame), 1), dtype=float)
    design = np.concatenate([intercept, *matrices], axis=1) if matrices else intercept
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise ValueError(
            f"{context}: nuisance design is rank deficient "
            f"(rank={rank}, columns={design.shape[1]})."
        )
    return design


def residualize_complete(
    values: np.ndarray,
    design: np.ndarray,
    *,
    label: str,
    context: str,
) -> np.ndarray:
    """Residualize a complete vector and reject degenerate adjusted values."""
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError(f"{context}: {label} contains non-finite values.")
    if design.shape[0] != array.shape[0]:
        raise ValueError(f"{context}: nuisance design row count does not match {label}.")
    beta, *_ = np.linalg.lstsq(design, array, rcond=None)
    residuals = array - design @ beta
    if not np.isfinite(residuals).all() or np.isclose(np.std(residuals), 0.0):
        raise ValueError(f"{context}: {label} has no usable variation after covariate adjustment.")
    return residuals
