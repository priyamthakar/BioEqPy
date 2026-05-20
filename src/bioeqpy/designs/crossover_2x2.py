"""Helpers for standard 2x2 crossover bioequivalence designs."""

from __future__ import annotations

import pandas as pd

from bioeqpy.core.exceptions import DesignError
from bioeqpy.core.types import BEStudy


def validate_2x2_completeness(study: BEStudy) -> pd.DataFrame:
    """Return a complete 2x2 dataframe or raise for invalid period/treatment rows."""
    frame = study.frame.copy()
    counts = frame.groupby("subject").size()
    if (counts < 2).any():
        missing = counts[counts < 2].index.astype(str).tolist()
        raise DesignError(f"Subjects with incomplete 2x2 records: {', '.join(missing)}")
    duplicated = frame.duplicated(subset=["subject", "period"])
    if duplicated.any():
        raise DesignError("Duplicate subject-period rows are not allowed.")
    invalid = frame.groupby("subject")["treatment"].nunique()
    if (invalid != 2).any():
        bad = invalid[invalid != 2].index.astype(str).tolist()
        raise DesignError(f"Each 2x2 subject must receive both T and R. Invalid subjects: {', '.join(bad)}")
    return frame

