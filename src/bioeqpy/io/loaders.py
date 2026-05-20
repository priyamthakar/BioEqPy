"""Load BE datasets from tabular files."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from bioeqpy.core.constants import REQUIRED_COLUMNS
from bioeqpy.core.types import BEStudy
from bioeqpy.designs import detect_design
from bioeqpy.io.validators import validate_dataset, validate_parameter


def load_dataset(path: str | Path) -> pd.DataFrame:
    """Load a CSV or Excel BE dataset."""
    source = Path(path)
    if source.suffix.lower() in {".xlsx", ".xls"}:
        data = pd.read_excel(source)
    else:
        data = pd.read_csv(source)
    data = data.rename(columns={column: str(column).strip() for column in data.columns})
    validate_dataset(data)
    return data


def infer_parameters(data: pd.DataFrame) -> list[str]:
    """Infer PK parameter columns from a BE dataset."""
    excluded = set(REQUIRED_COLUMNS)
    parameters: list[str] = []
    for column in data.columns:
        if column in excluded:
            continue
        values = pd.to_numeric(data[column], errors="coerce")
        if values.notna().all():
            parameters.append(column)
    return parameters


def load_study(path_or_data: str | Path | pd.DataFrame, parameter: str) -> BEStudy:
    """Load one PK parameter as a BEStudy with natural-log transformed values."""
    data = load_dataset(path_or_data) if not isinstance(path_or_data, pd.DataFrame) else path_or_data.copy()
    validate_dataset(data)
    validate_parameter(data, parameter)
    design = detect_design(data)
    raw_values = pd.to_numeric(data[parameter], errors="raise").to_numpy(dtype=float)
    return BEStudy(
        subjects=data["subject"].to_numpy(),
        sequences=data["sequence"].astype(str).str.upper().to_numpy(),
        periods=data["period"].astype(int).to_numpy(),
        treatments=data["treatment"].astype(str).str.upper().to_numpy(),
        values=np.log(raw_values),
        parameter_name=parameter,
        design=design,
        raw_values=raw_values,
        metadata={"source_rows": len(data)},
    )

