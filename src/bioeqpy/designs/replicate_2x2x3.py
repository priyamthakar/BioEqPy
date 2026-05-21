"""2x2x3 partial replicate crossover design (TRR/RTR/RRT)."""
from __future__ import annotations

import pandas as pd

from bioeqpy.core.types import BEStudy, DesignSpec

SPEC_2x2x3 = DesignSpec(
    name="2x2x3",
    n_periods=3,
    n_sequences=3,
    n_treatments=2,
    sequences=["RRT", "RTR", "TRR"],
    is_replicate=True,
    allows_scaled_abe=True,
)


def validate_2x2x3_completeness(study: BEStudy) -> pd.DataFrame:
    """Validate and return a complete 2x2x3 dataframe."""
    from bioeqpy.core.exceptions import DesignError

    frame = study.frame.copy()

    counts = frame.groupby("subject").size()
    if (counts < 3).any():
        missing = counts[counts < 3].index.astype(str).tolist()
        raise DesignError(f"Subjects with incomplete 2x2x3 records: {', '.join(missing)}")

    duplicated = frame.duplicated(subset=["subject", "period"])
    if duplicated.any():
        raise DesignError("Duplicate subject-period rows are not allowed.")

    # Each subject must receive T exactly once and R exactly twice
    trt_counts = frame.groupby(["subject", "treatment"]).size().unstack(fill_value=0)
    bad = ~((trt_counts.get("T", pd.Series(0, index=trt_counts.index)) == 1) &
            (trt_counts.get("R", pd.Series(0, index=trt_counts.index)) == 2))
    if bad.any():
        raise DesignError(
            "Each subject must receive T once and R twice in a 2x2x3 partial replicate design."
        )

    return frame
