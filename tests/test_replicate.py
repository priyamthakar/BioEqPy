"""Tests for the 2x2x4 full replicate crossover design."""
from pathlib import Path

from bioeqpy.io.loaders import load_study
from bioeqpy.anova import fit_anova
from bioeqpy import analyze

DATA = Path(__file__).parent / "reference_data" / "example_2x2x4.csv"


def test_2x2x4_design_detected():
    study = load_study(DATA, "AUClast")
    assert study.design.name == "2x2x4"
    assert study.design.is_replicate is True


def test_2x2x4_within_subject_variances():
    study = load_study(DATA, "AUClast")
    result = fit_anova(study)
    assert result.swr is not None
    assert result.swt is not None
    assert result.swr > 0
    assert result.swt > 0


def test_2x2x4_ci_is_reasonable():
    results = analyze(DATA, parameters=["AUClast"])
    r = results[0]
    assert r.design.name == "2x2x4"
    # The test values are ~7% higher than reference values, so CI should be around 107-108%
    assert 100 < r.ci.point_estimate < 115
    assert r.ci.lower > 70
    assert r.ci.upper < 140
