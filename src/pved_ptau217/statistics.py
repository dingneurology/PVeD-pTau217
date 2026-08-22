"""Statistical models with explicit safeguards and reproducible inference."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from .validation import require_columns, require_finite_variation, require_positive


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    standard_deviation = values.std(skipna=True, ddof=0)
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        return pd.Series(np.nan, index=series.index, dtype=float)
    return (values - values.mean(skipna=True)) / standard_deviation


def bh_fdr(series: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=series.index, dtype=float)
    valid = series.notna() & np.isfinite(series)
    if valid.any():
        adjusted.loc[valid] = multipletests(
            series.loc[valid].astype(float),
            method="fdr_bh",
        )[1]
    return adjusted


def _complete_case_design(
    frame: pd.DataFrame,
    *,
    outcome: str,
    ptau: str,
    pved: str,
    covariates: Sequence[str],
    categorical: Iterable[str],
) -> tuple[pd.Series, pd.DataFrame, dict[str, int]]:
    if outcome in covariates:
        raise ValueError(
            f"Invalid model: outcome '{outcome}' is also listed as a covariate."
        )
    required = [outcome, ptau, pved, *covariates]
    require_columns(frame, required, "model table")
    require_positive(frame, ptau, "model table", allow_missing=True)

    raw = frame[required].copy()
    missing_by_variable = {
        column: int(raw[column].isna().sum()) for column in required
    }
    complete = raw.dropna(subset=required).copy()
    complete[ptau] = np.log10(pd.to_numeric(complete[ptau], errors="raise"))
    numeric_for_check = [outcome, ptau, pved] + [
        column for column in covariates if column not in set(categorical)
    ]
    require_finite_variation(complete, numeric_for_check)

    response = zscore(complete[outcome]).rename("outcome_z")
    design = pd.DataFrame(index=complete.index)
    design["ptau_z"] = zscore(complete[ptau])
    design["pved_z"] = zscore(complete[pved])
    design["ptau_x_pved"] = design["ptau_z"] * design["pved_z"]

    categorical_set = set(categorical)
    for covariate in covariates:
        if covariate in categorical_set:
            encoded = pd.get_dummies(
                complete[covariate].astype("category"),
                prefix=covariate,
                drop_first=True,
                dtype=float,
            )
            if encoded.empty:
                raise ValueError(
                    f"Categorical covariate '{covariate}' has fewer than two levels."
                )
            design = pd.concat([design, encoded], axis=1)
        else:
            design[covariate] = zscore(complete[covariate])

    design = sm.add_constant(design, has_constant="add")
    if response.isna().any() or design.isna().any().any():
        raise ValueError("Unexpected missing value after complete-case preparation.")
    return response, design.astype(float), missing_by_variable


def fit_interaction_ols(
    frame: pd.DataFrame,
    *,
    outcome: str,
    ptau: str,
    pved: str,
    covariates: Sequence[str],
    categorical: Iterable[str],
    minimum_n: int = 30,
    covariance: str = "HC3",
    max_condition_number: float = 1_000_000,
) -> dict[str, object]:
    response, design, missing = _complete_case_design(
        frame,
        outcome=outcome,
        ptau=ptau,
        pved=pved,
        covariates=covariates,
        categorical=categorical,
    )
    n = len(response)
    if n < minimum_n:
        raise ValueError(f"Complete-case sample is {n}, below minimum {minimum_n}.")

    rank = int(np.linalg.matrix_rank(design.to_numpy()))
    if rank != design.shape[1]:
        raise ValueError(
            f"Design matrix is rank deficient: rank={rank}, columns={design.shape[1]}."
        )
    condition_number = float(np.linalg.cond(design.to_numpy()))
    if not np.isfinite(condition_number) or condition_number > max_condition_number:
        raise ValueError(
            f"Design condition number {condition_number:.3g} exceeds "
            f"{max_condition_number:.3g}."
        )

    fitted = sm.OLS(response, design).fit(cov_type=covariance)
    term = "ptau_x_pved"
    confidence = fitted.conf_int().loc[term]
    return {
        "n": n,
        "interaction_beta": float(fitted.params[term]),
        "robust_se": float(fitted.bse[term]),
        "ci_low": float(confidence.iloc[0]),
        "ci_high": float(confidence.iloc[1]),
        "p_value": float(fitted.pvalues[term]),
        "r_squared": float(fitted.rsquared),
        "adjusted_r_squared": float(fitted.rsquared_adj),
        "covariance": covariance,
        "matrix_rank": rank,
        "condition_number": condition_number,
        "missing_by_variable": missing,
    }


def run_outcome_family(
    frame: pd.DataFrame,
    *,
    cohort: str,
    marker_role: str,
    ptau: str,
    pved: str,
    outcomes: dict[str, dict[str, str]],
    model_sets: dict[str, Sequence[str]],
    categorical: Iterable[str],
    minimum_n: int,
    covariance: str,
    max_condition_number: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for outcome, metadata in outcomes.items():
        if outcome not in frame:
            continue
        for model_id, covariates in model_sets.items():
            if outcome in covariates:
                rows.append(
                    {
                        "cohort": cohort,
                        "marker_role": marker_role,
                        "outcome": outcome,
                        "outcome_label": metadata["label"],
                        "analysis_family": metadata["family"],
                        "model_id": model_id,
                        "model_status": "skipped_outcome_is_covariate",
                    }
                )
                continue
            try:
                result = fit_interaction_ols(
                    frame,
                    outcome=outcome,
                    ptau=ptau,
                    pved=pved,
                    covariates=list(covariates),
                    categorical=categorical,
                    minimum_n=minimum_n,
                    covariance=covariance,
                    max_condition_number=max_condition_number,
                )
                rows.append(
                    {
                        "cohort": cohort,
                        "marker_role": marker_role,
                        "outcome": outcome,
                        "outcome_label": metadata["label"],
                        "analysis_family": metadata["family"],
                        "model_id": model_id,
                        "covariates": ";".join(covariates),
                        "model_status": "ok",
                        **{k: v for k, v in result.items() if k != "missing_by_variable"},
                    }
                )
            except ValueError as error:
                rows.append(
                    {
                        "cohort": cohort,
                        "marker_role": marker_role,
                        "outcome": outcome,
                        "outcome_label": metadata["label"],
                        "analysis_family": metadata["family"],
                        "model_id": model_id,
                        "covariates": ";".join(covariates),
                        "model_status": "failed_validation",
                        "note": str(error),
                    }
                )

    results = pd.DataFrame(rows)
    if results.empty:
        return results
    results["fdr"] = np.nan
    valid = results["model_status"].eq("ok")
    grouping = ["cohort", "marker_role", "model_id", "analysis_family"]
    results.loc[valid, "fdr"] = (
        results.loc[valid]
        .groupby(grouping, dropna=False)["p_value"]
        .transform(bh_fdr)
    )
    return results
