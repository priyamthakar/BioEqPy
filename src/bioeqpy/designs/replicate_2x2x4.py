"""2x2x4 full replicate crossover design."""
from __future__ import annotations

import pandas as pd

from bioeqpy.core.types import BEStudy, DesignSpec

SPEC_2x2x4 = DesignSpec(
    name="2x2x4",
    n_periods=4,
    n_sequences=2,
    n_treatments=2,
    sequences=["RTTR", "TRRT"],
    is_replicate=True,
    allows_scaled_abe=True,
)


def validate_2x2x4_completeness(study: BEStudy) -> pd.DataFrame:
    """Validate and return a complete 2x2x4 dataframe."""
    from bioeqpy.core.exceptions import DesignError

    frame = study.frame.copy()

    counts = frame.groupby("subject").size()
    if (counts < 4).any():
        missing = counts[counts < 4].index.astype(str).tolist()
        raise DesignError(f"Subjects with incomplete 2x2x4 records: {', '.join(missing)}")

    duplicated = frame.duplicated(subset=["subject", "period"])
    if duplicated.any():
        raise DesignError("Duplicate subject-period rows are not allowed.")

    # Each subject must receive T exactly twice and R exactly twice
    trt_counts = frame.groupby(["subject", "treatment"]).size().unstack(fill_value=0)
    if not all((trt_counts.get("T", 0) == 2) & (trt_counts.get("R", 0) == 2)):
        raise DesignError("Each subject must receive T twice and R twice in a 2x2x4 design.")

    return frame
