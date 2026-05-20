"""Input validation for BioEqPy datasets."""

from __future__ import annotations

import numpy as np
import pandas as pd

from bioeqpy.core.constants import REQUIRED_COLUMNS
from bioeqpy.core.exceptions import ValidationError


def validate_dataset(data: pd.DataFrame) -> None:
    """Validate the common BE input schema."""
    missing = [column for column in REQUIRED_COLUMNS if column not in data.columns]
    if missing:
        raise ValidationError(f"Missing required columns: {', '.join(missing)}")
    if data.empty:
        raise ValidationError("Dataset is empty.")
    if data[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValidationError("Required design columns cannot contain missing values.")


def validate_parameter(data: pd.DataFrame, parameter: str) -> None:
    """Validate a numeric PK parameter column."""
    if parameter not in data.columns:
        raise ValidationError(f"Parameter column not found: {parameter}")
    values = pd.to_numeric(data[parameter], errors="coerce")
    if values.isna().any():
        raise ValidationError(f"Parameter column contains non-numeric or missing values: {parameter}")
    if np.any(values.to_numpy(dtype=float) <= 0):
        raise ValidationError(f"Parameter values must be positive before log transformation: {parameter}")

