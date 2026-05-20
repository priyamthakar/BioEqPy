"""Detect supported bioequivalence study designs."""

from __future__ import annotations

import pandas as pd

from bioeqpy.core.exceptions import DesignError
from bioeqpy.core.types import DesignSpec


def detect_design(data: pd.DataFrame) -> DesignSpec:
    """Detect the BE design from period, sequence, and treatment structure."""
    periods = sorted(data["period"].dropna().astype(int).unique().tolist())
    sequences = sorted(data["sequence"].dropna().astype(str).unique().tolist())
    treatments = sorted(data["treatment"].dropna().astype(str).unique().tolist())

    n_periods = len(periods)
    n_sequences = len(sequences)
    n_treatments = len(treatments)

    if n_periods == 2 and n_sequences == 2 and n_treatments == 2:
        expected = {tuple("TR"), tuple("RT")}
        observed = {tuple(str(seq).upper()) for seq in sequences}
        if observed != expected:
            raise DesignError("2x2 crossover requires TR and RT sequence labels.")
        return DesignSpec(
            name="2x2",
            n_periods=2,
            n_sequences=2,
            n_treatments=2,
            sequences=sequences,
            is_replicate=False,
            allows_scaled_abe=False,
        )

    if n_periods == 4 and n_sequences == 2 and n_treatments == 2:
        expected = {"RTTR", "TRRT"}
        observed = {str(seq).upper() for seq in sequences}
        if observed != expected:
            raise DesignError("2x2x4 full replicate crossover requires RTTR and TRRT sequence labels.")
        return DesignSpec(
            name="2x2x4",
            n_periods=4,
            n_sequences=2,
            n_treatments=2,
            sequences=sequences,
            is_replicate=True,
            allows_scaled_abe=True,
        )

    if n_periods == 1 and n_treatments == 2:
        return DesignSpec(
            name="parallel",
            n_periods=1,
            n_sequences=max(n_sequences, 1),
            n_treatments=2,
            sequences=sequences,
            is_replicate=False,
            allows_scaled_abe=False,
        )

    raise DesignError(
        "Unsupported design for v0.1. Supported: standard 2x2 crossover and basic parallel detection."
    )

