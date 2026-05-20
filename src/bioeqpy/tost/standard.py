"""Standard average bioequivalence TOST calculations."""

from __future__ import annotations

import math

from scipy import stats

from bioeqpy.core.constants import STANDARD_ACCEPTANCE
from bioeqpy.core.types import ANOVAResult, CIResult


def compute_ci(
    anova_result: ANOVAResult,
    alpha: float = 0.05,
    acceptance: tuple[float, float] = STANDARD_ACCEPTANCE,
    method: str = "standard",
) -> CIResult:
    """Compute the 100*(1-2*alpha)% CI on the geometric mean ratio."""
    tcrit = stats.t.ppf(1.0 - alpha, anova_result.residual_df)
    lower_log = anova_result.treatment_diff - tcrit * anova_result.se_diff
    upper_log = anova_result.treatment_diff + tcrit * anova_result.se_diff
    point = math.exp(anova_result.treatment_diff) * 100.0
    lower = math.exp(lower_log) * 100.0
    upper = math.exp(upper_log) * 100.0
    passed = acceptance[0] <= lower and upper <= acceptance[1]
    return CIResult(
        point_estimate=point,
        lower=lower,
        upper=upper,
        acceptance_lower=acceptance[0],
        acceptance_upper=acceptance[1],
        passed=passed,
        method=method,
    )

