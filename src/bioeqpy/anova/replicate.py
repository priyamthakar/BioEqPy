"""Within-subject variance component estimation for replicate designs."""
from __future__ import annotations

import numpy as np
import pandas as pd

from bioeqpy.core.types import BEStudy


def estimate_within_subject_variances(study: BEStudy) -> tuple[float, float, float, float]:
    """Estimate within-subject variances for T and R in a 2x2x4 design.

    Returns: (sWR, sWT, MSW_R, MSW_T)
    sWR and sWT are within-subject standard deviations.
    MSW_R and MSW_T are within-subject mean squares.
    """
    frame = study.frame.copy()

    results_r: list[float] = []
    results_t: list[float] = []

    for subject, subj_frame in frame.groupby("subject"):
        r_vals = subj_frame[subj_frame["treatment"] == "R"]["value"].to_numpy()
        t_vals = subj_frame[subj_frame["treatment"] == "T"]["value"].to_numpy()

        if len(r_vals) == 2:
            results_r.append((r_vals[0] - r_vals[1]) ** 2 / 2)
        if len(t_vals) == 2:
            results_t.append((t_vals[0] - t_vals[1]) ** 2 / 2)

    MSW_R = float(np.mean(results_r)) if results_r else 0.0
    MSW_T = float(np.mean(results_t)) if results_t else 0.0
    sWR = float(np.sqrt(MSW_R))
    sWT = float(np.sqrt(MSW_T))

    return sWR, sWT, MSW_R, MSW_T
