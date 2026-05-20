"""Tests for the 2x2x3 partial replicate crossover design (TRR/RTR/RRT)."""
from pathlib import Path

import numpy as np
import pytest

from bioeqpy import analyze
from bioeqpy.anova import fit_anova
from bioeqpy.io.loaders import load_study

DATA = Path(__file__).parent / "reference_data" / "example_2x2x3.csv"


@pytest.fixture(scope="module")
def study():
    return load_study(DATA, "AUClast")


@pytest.fixture(scope="module")
def anova_result(study):
    return fit_anova(study)


@pytest.fixture(scope="module")
def be_result():
    return analyze(DATA, parameters=["AUClast"])[0]


# --- Design detection ---

def test_2x2x3_design_detected(study):
    assert study.design.name == "2x2x3"


def test_2x2x3_design_spec(study):
    d = study.design
    assert d.n_periods == 3
    assert d.n_sequences == 3
    assert d.n_treatments == 2
    assert d.is_replicate is True
    assert d.allows_scaled_abe is True
    assert set(d.sequences) == {"TRR", "RTR", "RRT"}


# --- Within-subject variance ---

def test_2x2x3_swr_is_estimable(anova_result):
    """sWR must be positive (R appears twice per subject)."""
    assert anova_result.swr is not None
    assert anova_result.swr > 0


def test_2x2x3_swt_is_none(anova_result):
    """sWT is not estimable from a partial replicate — must be None."""
    assert anova_result.swt is None


# --- ANOVA structure ---

def test_2x2x3_residual_df(study, anova_result):
    """Residual df must equal 2n - 3 for a 3-period partial replicate."""
    n = int(study.frame["subject"].nunique())
    expected_df = 2 * n - 3
    assert anova_result.residual_df == expected_df


def test_2x2x3_source_table_sources(anova_result):
    sources = set(anova_result.source_table["Source"])
    required = {"Sequence", "Subject within Sequence", "Period", "Treatment", "Residual"}
    assert required.issubset(sources)


def test_2x2x3_source_table_no_nan_ss(anova_result):
    ss = anova_result.source_table["SS"]
    assert ss.notna().all(), f"NaN SS found:\n{anova_result.source_table}"


def test_2x2x3_source_table_df_sum(study, anova_result):
    """Total DF must equal 3n - 1."""
    n = int(study.frame["subject"].nunique())
    expected = 3 * n - 1
    total = int(anova_result.source_table["DF"].sum())
    assert total == expected, f"Expected {expected} total DF, got {total}"


def test_2x2x3_source_table_ss_partitions_total(study, anova_result):
    """SS must partition SS_total exactly (sequential decomposition)."""
    frame = study.frame
    ss_total = float(((frame["value"] - frame["value"].mean()) ** 2).sum())
    ss_sum = float(anova_result.source_table["SS"].sum())
    assert abs(ss_sum - ss_total) < 1e-8, (
        f"SS partition fails: table={ss_sum:.10f}, total={ss_total:.10f}, "
        f"delta={abs(ss_sum - ss_total):.2e}"
    )


def test_2x2x3_sequence_f_uses_subseq_denominator(anova_result):
    """Sequence F must use Subject(Sequence) MS as denominator (split-plot)."""
    tbl = anova_result.source_table.set_index("Source")
    ms_seq = float(tbl.loc["Sequence", "MS"])
    ms_sub_seq = float(tbl.loc["Subject within Sequence", "MS"])
    f_seq = float(tbl.loc["Sequence", "F"])
    expected = ms_seq / ms_sub_seq
    assert abs(f_seq - expected) < 1e-10


def test_2x2x3_residual_ms_matches(anova_result):
    tbl = anova_result.source_table.set_index("Source")
    assert abs(tbl.loc["Residual", "MS"] - anova_result.residual_ms) < 1e-12


# --- CI and BE decision ---

def test_2x2x3_ci_is_reasonable(be_result):
    """GMR should be around 104% and CI within 80–125% for these data."""
    assert 100.0 < be_result.ci.point_estimate < 115.0
    assert be_result.ci.lower > 80.0
    assert be_result.ci.upper < 125.0


def test_2x2x3_bioequivalent(be_result):
    assert be_result.ci.passed
    assert be_result.conclusion == "Bioequivalent"


def test_2x2x3_abel_method():
    """ABEL method uses sWR; with low sWR falls back to 80–125% limits."""
    results = analyze(DATA, parameters=["AUClast"], method="abel")
    r = results[0]
    assert r.ci.method == "abel"
    assert r.ci.acceptance_lower == pytest.approx(80.0)
    assert r.ci.acceptance_upper == pytest.approx(125.0)
    assert r.ci.passed


def test_2x2x3_rsabe_method():
    """RSABE method must run without error and return a CI."""
    results = analyze(DATA, parameters=["AUClast"], method="rsabe")
    r = results[0]
    assert r.ci.method == "rsabe"
    assert r.ci.point_estimate > 0
