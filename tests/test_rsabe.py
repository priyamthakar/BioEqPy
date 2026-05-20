from pathlib import Path
import math
import numpy as np
import pytest
from bioeqpy import analyze
from bioeqpy.anova import fit_anova
from bioeqpy.io.loaders import load_study
from bioeqpy.tost.rsabe import (
    compute_abel_ci,
    compute_rsabe_ci,
    SWR_THRESHOLD,
    ABEL_K,
    ABEL_CAP_LOWER,
    ABEL_CAP_UPPER,
)

DATA_4 = Path(__file__).parent / "reference_data" / "example_2x2x4.csv"


def test_abel_method_label():
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    result = compute_abel_ci(anova)
    assert result.method == "abel"


def test_abel_limits_when_swr_below_threshold():
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    if anova.swr < SWR_THRESHOLD:
        result = compute_abel_ci(anova)
        assert result.acceptance_lower == 80.0
        assert result.acceptance_upper == 125.0


def test_abel_limits_when_swr_above_threshold():
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    if anova.swr >= SWR_THRESHOLD:
        result = compute_abel_ci(anova)
        expected_hi = min(math.exp(ABEL_K * anova.swr) * 100.0, ABEL_CAP_UPPER)
        assert abs(result.acceptance_upper - expected_hi) < 0.01


def test_abel_limits_always_within_cap():
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    result = compute_abel_ci(anova)
    assert result.acceptance_lower >= ABEL_CAP_LOWER
    assert result.acceptance_upper <= ABEL_CAP_UPPER


def test_rsabe_method_label():
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    result = compute_rsabe_ci(anova)
    assert result.method == "rsabe"


def test_rsabe_passed_is_bool():
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    result = compute_rsabe_ci(anova)
    assert isinstance(result.passed, bool)


def test_analyze_abel():
    results = analyze(DATA_4, parameters=["AUClast"], method="abel")
    assert results[0].ci.method == "abel"
    assert results[0].conclusion in {"Bioequivalent", "Not bioequivalent"}


def test_analyze_rsabe():
    results = analyze(DATA_4, parameters=["AUClast"], method="rsabe")
    assert results[0].ci.method == "rsabe"
    assert results[0].conclusion in {"Bioequivalent", "Not bioequivalent"}


def test_se_diff_corrected():
    """SE should be sqrt(pooled_MSW / n), not sqrt(2 * pooled_MSW / n)."""
    from bioeqpy.anova.replicate import estimate_within_subject_variances
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    _, _, MSW_R, MSW_T = estimate_within_subject_variances(study)
    pooled_MSW = (MSW_R + MSW_T) / 2.0
    n = study.frame["subject"].nunique()
    expected_se = math.sqrt(pooled_MSW / n)
    assert abs(anova.se_diff - expected_se) < 1e-10


def test_2x2x4_source_table_no_nan_treatment():
    """Treatment row in 2x2x4 source table must not contain NaN values."""
    study = load_study(DATA_4, "AUClast")
    anova = fit_anova(study)
    trt_row = anova.source_table[anova.source_table["Source"] == "Treatment"].iloc[0]
    assert not np.isnan(trt_row["SS"])
    assert not np.isnan(trt_row["MS"])
    assert not np.isnan(trt_row["F"])
    assert not np.isnan(trt_row["p-value"])


def test_rsabe_regulatory_basis():
    results = analyze(DATA_4, parameters=["AUClast"], method="rsabe")
    assert "RSABE" in results[0].regulatory_basis or "rsabe" in results[0].regulatory_basis.lower()


def test_abel_regulatory_basis():
    results = analyze(DATA_4, parameters=["AUClast"], method="abel")
    assert "ABEL" in results[0].regulatory_basis or "abel" in results[0].regulatory_basis.lower()
