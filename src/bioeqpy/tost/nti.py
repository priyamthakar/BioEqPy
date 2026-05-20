"""NTI-specific bioequivalence CI computation."""
from __future__ import annotations

from bioeqpy.core.constants import NTI_ACCEPTANCE
from bioeqpy.core.types import ANOVAResult, CIResult


def compute_nti_ci(
    anova_result: ANOVAResult,
    alpha: float = 0.05,
    acceptance: tuple[float, float] = NTI_ACCEPTANCE,
) -> CIResult:
    """Compute 90% CI with NTI tightened acceptance limits.

    EMA tightened limits: 90.00-111.11%.
    Some product-specific FDA guidances use 80-125% with an additional
    intra-subject variability ratio constraint -- this function implements
    the EMA tightened-limits approach only.
    """
    from bioeqpy.tost.standard import compute_ci

    return compute_ci(anova_result, alpha=alpha, acceptance=acceptance, method="nti")
