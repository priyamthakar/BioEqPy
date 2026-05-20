from pathlib import Path
import pytest
from bioeqpy import analyze
from bioeqpy.io.loaders import load_study
from bioeqpy.diagnostics.carryover import carryover_effect

DATA = Path(__file__).parent / "reference_data" / "example_2x2.csv"


def test_nti_tighter_limits_fail_standard_study():
    """The example dataset passes 80-125% but may fail 90-111.11%."""
    results = analyze(DATA, parameters=["AUClast"], method="nti")
    r = results[0]
    # NTI limits are tighter -- check that limits are set correctly
    assert r.ci.acceptance_lower == 90.0
    assert r.ci.acceptance_upper == 111.11
    # The method should be labeled nti
    assert r.ci.method == "nti"


def test_nti_conclusion_label():
    """NTI analysis should return a valid conclusion string."""
    results = analyze(DATA, parameters=["AUClast"], method="nti")
    assert results[0].conclusion in {"Bioequivalent", "Not bioequivalent"}


def test_carryover_returns_pvalue():
    study = load_study(DATA, "AUClast")
    result = carryover_effect(study.frame)
    assert "p_value" in result
    assert result["p_value"] is not None
    assert 0.0 <= result["p_value"] <= 1.0


def test_carryover_in_diagnostics():
    """check_assumptions should now include carryover results."""
    results = analyze(DATA, parameters=["AUClast"])
    diag = results[0].diagnostics
    assert "carryover" in diag
    assert "p_value" in diag["carryover"]


def test_standard_analysis_unaffected():
    """Existing standard analysis must still pass all assertions."""
    results = analyze(DATA, parameters=["AUClast"])
    r = results[0]
    assert r.ci.passed
    assert r.ci.acceptance_lower == 80.0
    assert r.ci.acceptance_upper == 125.0
