"""Typed data containers for bioequivalence analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class DesignSpec:
    """Crossover or parallel design specification."""

    name: str
    n_periods: int
    n_sequences: int
    n_treatments: int
    sequences: list[str]
    is_replicate: bool = False
    allows_scaled_abe: bool = False


@dataclass
class BEStudy:
    """Container for one log-transformed PK parameter dataset."""

    subjects: np.ndarray
    sequences: np.ndarray
    periods: np.ndarray
    treatments: np.ndarray
    values: np.ndarray
    parameter_name: str
    design: DesignSpec
    raw_values: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def frame(self) -> pd.DataFrame:
        """Return the study as an analysis dataframe."""
        return pd.DataFrame(
            {
                "subject": self.subjects,
                "sequence": self.sequences,
                "period": self.periods,
                "treatment": self.treatments,
                "value": self.values,
                "raw_value": self.raw_values if self.raw_values is not None else np.exp(self.values),
            }
        )


@dataclass
class ANOVAResult:
    """ANOVA decomposition output for one PK parameter."""

    source_table: pd.DataFrame
    residual_ms: float
    residual_df: float
    treatment_diff: float
    se_diff: float
    residuals: pd.Series
    fitted: pd.Series
    swr: float | None = None
    swt: float | None = None


@dataclass(frozen=True)
class CIResult:
    """Confidence interval on the geometric mean ratio."""

    point_estimate: float
    lower: float
    upper: float
    acceptance_lower: float
    acceptance_upper: float
    passed: bool
    method: str


@dataclass
class BEResult:
    """Complete bioequivalence analysis result for one PK parameter."""

    parameter_name: str
    n_subjects: int
    design: DesignSpec
    summary_stats: pd.DataFrame
    anova: ANOVAResult
    ci: CIResult
    diagnostics: dict[str, Any]
    conclusion: str
    regulatory_basis: str

