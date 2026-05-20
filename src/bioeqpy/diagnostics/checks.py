"""Diagnostic checks for BE analyses."""

from __future__ import annotations

import numpy as np
from scipy import stats

from bioeqpy.core.types import ANOVAResult, BEStudy


def check_assumptions(study: BEStudy, anova_result: ANOVAResult) -> dict[str, object]:
    """Run lightweight diagnostics for the fitted model."""
    residuals = anova_result.residuals.dropna().to_numpy(dtype=float)
    diagnostics: dict[str, object] = {
        "n_residuals": int(len(residuals)),
        "warnings": [],
    }
    if len(residuals) >= 3:
        shapiro = stats.shapiro(residuals)
        diagnostics["normality_p_value"] = float(shapiro.pvalue)
        if shapiro.pvalue < 0.05:
            diagnostics["warnings"].append("Residual normality test is significant at p < 0.05.")
    else:
        diagnostics["normality_p_value"] = None
        diagnostics["warnings"].append("Too few residuals for Shapiro-Wilk normality testing.")

    std = float(np.std(residuals, ddof=1)) if len(residuals) > 1 else 0.0
    if std > 0:
        studentized = residuals / std
        diagnostics["outlier_indices"] = np.where(np.abs(studentized) > 3.0)[0].tolist()
    else:
        diagnostics["outlier_indices"] = []
    diagnostics["design"] = study.design.name

    if study.design.name == "2x2":
        from bioeqpy.diagnostics.carryover import carryover_effect

        carryover = carryover_effect(study.frame)
        diagnostics["carryover"] = carryover
        if carryover.get("significant"):
            diagnostics["warnings"].append(
                f"Carryover test significant: p = {carryover['p_value']:.4f}. "
                "Interpret with caution."
            )

    return diagnostics

