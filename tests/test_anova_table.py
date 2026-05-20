"""Structural validation tests for the 2x2 ANOVA source table."""

from pathlib import Path

import numpy as np
import pytest

from bioeqpy import analyze
from bioeqpy.anova import fit_anova
from bioeqpy.io.loaders import load_study

DATA = Path(__file__).parent / "reference_data" / "example_2x2.csv"


@pytest.fixture(scope="module")
def auclast_result():
    return analyze(DATA, parameters=["AUClast"])[0]


@pytest.fixture(scope="module")
def auclast_anova():
    study = load_study(DATA, "AUClast")
    return fit_anova(study), study


def test_source_table_has_no_nan_ss(auclast_result):
    """All SS values must be finite — no NaN placeholders."""
    ss_col = auclast_result.anova.source_table["SS"]
    assert ss_col.notna().all(), f"NaN SS found:\n{auclast_result.anova.source_table}"


def test_source_table_df_sums_to_n_minus_1(auclast_result):
    """Total DF = n_observations - 1 = 16 - 1 = 15."""
    total_df = int(auclast_result.anova.source_table["DF"].sum())
    assert total_df == 15, f"Expected 15 total DF, got {total_df}"


def test_source_table_ss_partitions_total(auclast_anova):
    """SS_seq + SS_sub_seq + SS_period + SS_treatment + SS_residual == SS_total."""
    anova_result, study = auclast_anova
    frame = study.frame
    ss_total = float(((frame["value"] - frame["value"].mean()) ** 2).sum())
    ss_sum = float(anova_result.source_table["SS"].sum())
    assert abs(ss_sum - ss_total) < 1e-8, (
        f"SS partition fails: table sum={ss_sum:.10f}, total={ss_total:.10f}, "
        f"delta={abs(ss_sum - ss_total):.2e}"
    )


def test_sequence_f_uses_subject_sequence_denominator(auclast_result):
    """Sequence F must equal MS_seq / MS_sub_seq, not MS_seq / MS_residual."""
    tbl = auclast_result.anova.source_table.set_index("Source")
    ms_seq = float(tbl.loc["Sequence", "MS"])
    ms_sub_seq = float(tbl.loc["Subject within Sequence", "MS"])
    f_seq = float(tbl.loc["Sequence", "F"])
    expected = ms_seq / ms_sub_seq
    assert abs(f_seq - expected) < 1e-10, (
        f"Sequence F={f_seq:.6f} != MS_seq/MS_sub_seq={expected:.6f}. "
        "Check split-plot denominator."
    )


def test_residual_ms_matches_anova_result(auclast_result):
    """Residual MS in source table must match ANOVAResult.residual_ms."""
    tbl = auclast_result.anova.source_table.set_index("Source")
    table_ms = float(tbl.loc["Residual", "MS"])
    assert abs(table_ms - auclast_result.anova.residual_ms) < 1e-12, (
        f"table residual MS={table_ms} != anova.residual_ms={auclast_result.anova.residual_ms}"
    )


def test_ci_unchanged_after_anova_fix(auclast_result):
    """Treatment diff/SE (and hence CI) must be unaffected by the SS fix."""
    assert auclast_result.ci.passed
    assert 100.0 < auclast_result.ci.point_estimate < 104.0
    assert auclast_result.ci.lower > 80.0
    assert auclast_result.ci.upper < 125.0


def test_all_parameters_have_no_nan_ss():
    """Verify no NaN SS for all three PK parameters."""
    results = analyze(DATA)
    for r in results:
        ss_col = r.anova.source_table["SS"]
        assert ss_col.notna().all(), (
            f"NaN SS in {r.parameter_name}:\n{r.anova.source_table}"
        )
