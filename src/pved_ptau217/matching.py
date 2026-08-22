"""Deterministic closest-date matching for external validation."""

from __future__ import annotations

import pandas as pd

from .validation import require_columns, require_unique, require_valid_dates


def closest_row_within_window(
    index: pd.DataFrame,
    measurements: pd.DataFrame,
    *,
    participant: str,
    index_date: str,
    measurement_date: str,
    required_value: str,
    window_months: float,
    prefix: str,
    tie_break: str = "earlier",
) -> pd.DataFrame:
    """Select one eligible measurement row nearest each verified index date.

    The index table must contain exactly one row per participant. Measurement
    rows with missing required values or dates are ineligible. No earliest-date
    fallback is performed.
    """
    if tie_break not in {"earlier", "later"}:
        raise ValueError("tie_break must be 'earlier' or 'later'.")
    if window_months <= 0:
        raise ValueError("window_months must be positive.")

    require_columns(index, [participant, index_date], "index")
    require_columns(
        measurements,
        [participant, measurement_date, required_value],
        "measurements",
    )
    require_unique(index, participant, "index")
    left = require_valid_dates(index, [index_date], "index", allow_missing=False)
    right = require_valid_dates(
        measurements,
        [measurement_date],
        "measurements",
        allow_missing=True,
    )
    right = right.loc[
        right[measurement_date].notna() & right[required_value].notna()
    ].copy()

    candidates = left[[participant, index_date]].merge(
        right,
        on=participant,
        how="left",
        validate="one_to_many",
    )
    candidates["signed_months"] = (
        candidates[measurement_date] - candidates[index_date]
    ).dt.days / 30.4375
    candidates["absolute_months"] = candidates["signed_months"].abs()
    eligible = candidates.loc[
        candidates["absolute_months"].le(float(window_months))
    ].copy()

    if eligible.empty:
        return pd.DataFrame(columns=[participant])

    tie_sign = 1 if tie_break == "earlier" else -1
    eligible["_tie_order"] = eligible["signed_months"] * tie_sign
    eligible = eligible.sort_values(
        [participant, "absolute_months", "_tie_order", measurement_date],
        kind="mergesort",
    )
    selected = eligible.drop_duplicates(participant, keep="first").copy()

    protected = {participant, index_date, "_tie_order"}
    rename = {
        column: f"{prefix}_{column}"
        for column in selected.columns
        if column not in protected
    }
    selected = selected.drop(columns=["_tie_order"]).rename(columns=rename)
    return selected.drop(columns=[index_date], errors="ignore").reset_index(drop=True)
