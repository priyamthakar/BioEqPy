"""Numerical validation of TOST CI against reference-computed values.

These tests assert that the BioEqPy ANOVA and CI engine produces results
within 0.01 percentage points of values independently derived from
first-principles formulas applied to the same input data. They also lock
the current output for example_2x2.csv to catch future regressions.
"""
from pathlib import Path

import numpy as np
import pytest

from bioeqpy import analyze
from bioeqpy.anova import fit_anova
from bioeqpy.io.loaders import load_study

DATA = Path(__file__).parent / "reference_data" / "example_2x2.csv"
DATA_V = Path(__file__).parent / "reference_data" / "validation_grizzle.csv"

# ---------------------------------------------------------------------------
# Locked regression guards for example_2x2.csv
# Computed 2026-05-20 from the current implementation.
# If these fail, the core ANOVA or TOST logic has regressed.
# ---------------------------------------------------------------------------
LOCKED_LOWER = 100.5743676264
LOCKED_UPPER = 103.4850780291
LOCKED_POINT = 102.0193426833


def test_example_2x2_ci_bounds_locked():
    r = analyze(DATA, parameters=["AUClast"])[0]
    assert abs(r.ci.lower - LOCKED_LOWER) < 0.0001, (
        f"CI lower regressed: got {r.ci.lower:.6f}, expected {LOCKED_LOWER:.6f}"
    )
    assert abs(r.ci.upper - LOCKED_UPPER) < 0.0001, (
        f"CI upper regressed: got {r.ci.upper:.6f}, expected {LOCKED_UPPER:.6f}"
    )


def test_example_2x2_point_estimate_locked():
    r = analyze(DATA, parameters=["AUClast"])[0]
    assert abs(r.ci.point_estimate - LOCKED_POINT) < 0.0001, (
        f"GMR regressed: got {r.ci.point_estimate:.6f}, expected {LOCKED_POINT:.6f}"
    )


# ---------------------------------------------------------------------------
# First-principles validation against Grizzle-style reference dataset
# 20 subjects, 2x2 crossover, 10 per sequence (TR/RT).
# Reference values computed independently using the contrast formula:
#   D_i = beta0 + sign_i * beta1 + eps,  where sign = +1 for TR, -1 for RT
#   SE = sqrt(cov_beta[1,1]),  residual_df = n - 2
#   90% CI: exp(beta1 ± t(0.95, df) * SE) * 100
# ---------------------------------------------------------------------------
REF_TREATMENT_DIFF = 0.067183229225
REF_SE_DIFF = 0.003786061281
REF_RESIDUAL_DF = 18
REF_LOWER = 106.2492919872
REF_UPPER = 107.6536022965
REF_POINT = 106.9491422307


def test_grizzle_treatment_diff():
    study = load_study(DATA_V, "AUClast")
    anova = fit_anova(study)
    assert abs(anova.treatment_diff - REF_TREATMENT_DIFF) < 1e-8, (
        f"treatment_diff mismatch: {anova.treatment_diff:.12f} vs {REF_TREATMENT_DIFF:.12f}"
    )


def test_grizzle_se_diff():
    study = load_study(DATA_V, "AUClast")
    anova = fit_anova(study)
    assert abs(anova.se_diff - REF_SE_DIFF) < 1e-8, (
        f"se_diff mismatch: {anova.se_diff:.12f} vs {REF_SE_DIFF:.12f}"
    )


def test_grizzle_residual_df():
    study = load_study(DATA_V, "AUClast")
    anova = fit_anova(study)
    assert anova.residual_df == REF_RESIDUAL_DF


def test_grizzle_ci_lower():
    r = analyze(DATA_V, parameters=["AUClast"])[0]
    assert abs(r.ci.lower - REF_LOWER) < 0.01, (
        f"90%% CI lower: got {r.ci.lower:.6f}, reference {REF_LOWER:.6f}"
    )


def test_grizzle_ci_upper():
    r = analyze(DATA_V, parameters=["AUClast"])[0]
    assert abs(r.ci.upper - REF_UPPER) < 0.01, (
        f"90%% CI upper: got {r.ci.upper:.6f}, reference {REF_UPPER:.6f}"
    )


def test_grizzle_point_estimate():
    r = analyze(DATA_V, parameters=["AUClast"])[0]
    assert abs(r.ci.point_estimate - REF_POINT) < 0.01, (
        f"GMR: got {r.ci.point_estimate:.6f}, reference {REF_POINT:.6f}"
    )


def test_grizzle_passes_be():
    r = analyze(DATA_V, parameters=["AUClast"])[0]
    assert r.ci.passed
    assert r.conclusion == "Bioequivalent"
