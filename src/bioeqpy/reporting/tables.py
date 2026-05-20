"""Regulatory-style table builders."""

from __future__ import annotations

import numpy as np
import pandas as pd

from bioeqpy.core.types import BEStudy, BEResult


def summary_statistics(study: BEStudy) -> pd.DataFrame:
    """Compute per-treatment descriptive statistics on raw and log scales."""
    frame = study.frame
    rows = []
    for treatment, group in frame.groupby("treatment", sort=True):
        raw = group["raw_value"].to_numpy(dtype=float)
        log_values = group["value"].to_numpy(dtype=float)
        rows.append(
            {
                "Treatment": treatment,
                "N": int(len(group)),
                "Arithmetic Mean": float(np.mean(raw)),
                "SD": float(np.std(raw, ddof=1)) if len(raw) > 1 else 0.0,
                "CV%": float(np.std(raw, ddof=1) / np.mean(raw) * 100.0) if len(raw) > 1 else 0.0,
                "Geometric Mean": float(np.exp(np.mean(log_values))),
            }
        )
    return pd.DataFrame(rows)


def format_ci_table(results: list[BEResult]) -> pd.DataFrame:
    """Return a compact BE assessment table for multiple parameters."""
    return pd.DataFrame(
        [
            {
                "Parameter": result.parameter_name,
                "N": result.n_subjects,
                "GMR(%)": round(result.ci.point_estimate, 2),
                "90% CI": f"{result.ci.lower:.2f}-{result.ci.upper:.2f}",
                "Limits": f"{result.ci.acceptance_lower:.2f}-{result.ci.acceptance_upper:.2f}",
                "Conclusion": "BE" if result.ci.passed else "Not BE",
            }
            for result in results
        ]
    )

