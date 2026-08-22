"""Input validation that fails early on ambiguous analysis data."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd


def require_columns(frame: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    required = list(dict.fromkeys(columns))
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def require_unique(frame: pd.DataFrame, column: str, table_name: str) -> None:
    if frame[column].isna().any():
        raise ValueError(f"{table_name}.{column} contains missing values.")
    duplicated = frame[column].duplicated(keep=False)
    if duplicated.any():
        raise ValueError(
            f"{table_name}.{column} must be unique; "
            f"{int(duplicated.sum())} rows are duplicated."
        )


def require_valid_dates(
    frame: pd.DataFrame,
    columns: Iterable[str],
    table_name: str,
    *,
    allow_missing: bool,
) -> pd.DataFrame:
    checked = frame.copy()
    for column in columns:
        parsed = pd.to_datetime(checked[column], errors="coerce")
        invalid = checked[column].notna() & parsed.isna()
        if invalid.any():
            raise ValueError(
                f"{table_name}.{column} contains {int(invalid.sum())} invalid dates."
            )
        if not allow_missing and parsed.isna().any():
            raise ValueError(f"{table_name}.{column} contains missing dates.")
        checked[column] = parsed
    return checked


def require_positive(
    frame: pd.DataFrame,
    column: str,
    table_name: str,
    *,
    allow_missing: bool = True,
) -> None:
    values = pd.to_numeric(frame[column], errors="coerce")
    invalid_text = frame[column].notna() & values.isna()
    if invalid_text.any():
        raise ValueError(f"{table_name}.{column} contains nonnumeric values.")
    if not allow_missing and values.isna().any():
        raise ValueError(f"{table_name}.{column} contains missing values.")
    nonpositive = values.notna() & values.le(0)
    if nonpositive.any():
        raise ValueError(
            f"{table_name}.{column} contains {int(nonpositive.sum())} "
            "non-positive values; log10 transformation is undefined."
        )


def require_finite_variation(frame: pd.DataFrame, columns: Iterable[str]) -> None:
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        finite = values[np.isfinite(values)]
        if finite.nunique() < 2:
            raise ValueError(f"{column} has fewer than two finite values.")
