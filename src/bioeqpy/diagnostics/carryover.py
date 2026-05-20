"""Carryover effect test for 2x2 crossover designs."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def carryover_effect(frame: pd.DataFrame) -> dict[str, object]:
    """Test for first-order carryover in a 2x2 crossover design.

    Uses the standard sequence-by-period interaction test:
    compares mean period differences between TR and RT sequence groups.
    A significant result suggests possible carryover from period 1 to period 2.

    Returns a dict with keys:
    - carryover_diff: float, difference in mean period differences between sequences
    - t_stat: float, t-statistic
    - p_value: float, two-sided p-value
    - df: int, degrees of freedom
    - significant: bool, True if p < 0.05
    - interpretation: str, plain-language note
    """
    pivot = frame.pivot(index="subject", columns="period", values="value")
    sequences = frame.groupby("subject")["sequence"].first().str.upper()

    period_diffs = pivot[1] - pivot[2]

    tr_diffs = period_diffs[sequences == "TR"].dropna()
    rt_diffs = period_diffs[sequences == "RT"].dropna()

    if len(tr_diffs) < 2 or len(rt_diffs) < 2:
        return {
            "carryover_diff": None,
            "t_stat": None,
            "p_value": None,
            "df": None,
            "significant": False,
            "interpretation": "Insufficient data for carryover test",
        }

    # Two-sample t-test on period differences between sequence groups
    t_stat, p_value = stats.ttest_ind(tr_diffs, rt_diffs, equal_var=False)
    df = len(tr_diffs) + len(rt_diffs) - 2
    carryover_diff = float(tr_diffs.mean() - rt_diffs.mean())
    significant = bool(p_value < 0.05)

    if significant:
        interpretation = "Possible carryover detected. Interpret treatment comparison with caution."
    else:
        interpretation = "No significant carryover detected at alpha 0.05."

    return {
        "carryover_diff": round(carryover_diff, 6),
        "t_stat": round(float(t_stat), 4),
        "p_value": round(float(p_value), 4),
        "df": df,
        "significant": significant,
        "interpretation": interpretation,
    }
