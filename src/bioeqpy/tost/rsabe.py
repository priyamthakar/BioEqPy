"""EMA ABEL and FDA RSABE scaled bioequivalence methods."""
from __future__ import annotations

import math

from scipy import stats

from bioeqpy.core.types import ANOVAResult, CIResult

SWR_THRESHOLD = 0.294
ABEL_K = math.log(1.25) / 0.25
ABEL_CAP_LOWER = 69.84
ABEL_CAP_UPPER = 143.19
RSABE_THETA = (math.log(1.25)) ** 2 / (0.25 ** 2)


def compute_abel_ci(anova_result: ANOVAResult, alpha: float = 0.05) -> CIResult:
    """EMA Average Bioequivalence with Expanding Limits.

    When sWR >= 0.294, acceptance limits expand to [exp(-k*sWR), exp(k*sWR)]*100
    capped at 69.84-143.19%. Falls back to 80-125 when sWR < threshold.
    EMA BE guideline 2010/2018.
    """
    sWR = anova_result.swr
    if sWR is None:
        raise ValueError("ABEL requires sWR from a replicate design.")

    if sWR < SWR_THRESHOLD:
        lo, hi = 80.0, 125.0
    else:
        lo = max(math.exp(-ABEL_K * sWR) * 100.0, ABEL_CAP_LOWER)
        hi = min(math.exp(ABEL_K * sWR) * 100.0, ABEL_CAP_UPPER)

    tcrit = stats.t.ppf(1.0 - alpha, anova_result.residual_df)
    point = math.exp(anova_result.treatment_diff) * 100.0
    lower = math.exp(anova_result.treatment_diff - tcrit * anova_result.se_diff) * 100.0
    upper = math.exp(anova_result.treatment_diff + tcrit * anova_result.se_diff) * 100.0
    passed = lo <= lower and upper <= hi

    return CIResult(
        point_estimate=point,
        lower=lower,
        upper=upper,
        acceptance_lower=lo,
        acceptance_upper=hi,
        passed=passed,
        method="abel",
    )


def compute_rsabe_ci(anova_result: ANOVAResult, alpha: float = 0.05) -> CIResult:
    """FDA Reference-Scaled ABE for highly variable drugs.

    Tests linearized criterion H = (mu_T - mu_R)^2 - theta * sigma_WR^2 <= 0
    using Howe (1974) approximation for the upper 95% confidence bound.
    FDA draft guidance 2011.
    """
    sWR = anova_result.swr
    if sWR is None:
        raise ValueError("RSABE requires sWR from a replicate design.")

    d = anova_result.treatment_diff
    se_d = anova_result.se_diff
    df = anova_result.residual_df
    s2_WR = sWR ** 2

    tcrit_90 = stats.t.ppf(0.95, df)
    point = math.exp(d) * 100.0
    lower = math.exp(d - tcrit_90 * se_d) * 100.0
    upper = math.exp(d + tcrit_90 * se_d) * 100.0

    if sWR < SWR_THRESHOLD:
        lo, hi = 80.0, 125.0
        passed = lo <= lower and upper <= hi
    else:
        t95 = stats.t.ppf(0.95, df)
        A1 = 4.0 * d ** 2 * (se_d * t95) ** 2
        A2 = 2.0 * (RSABE_THETA * s2_WR) ** 2 / df
        ucb_H = d ** 2 - RSABE_THETA * s2_WR + math.sqrt(A1 + A2)
        passed = ucb_H <= 0.0 and 80.0 <= lower and upper <= 125.0
        lo, hi = 80.0, 125.0

    return CIResult(
        point_estimate=point,
        lower=lower,
        upper=upper,
        acceptance_lower=lo,
        acceptance_upper=hi,
        passed=passed,
        method="rsabe",
    )
