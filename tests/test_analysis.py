from pathlib import Path

from bioeqpy import analyze


DATA = Path(__file__).parent / "reference_data" / "example_2x2.csv"


def test_analyze_returns_all_numeric_parameters():
    results = analyze(DATA)

    assert [result.parameter_name for result in results] == ["AUClast", "AUCinf", "Cmax"]
    assert all(result.design.name == "2x2" for result in results)
    assert all(result.anova.residual_df == 6 for result in results)


def test_auclast_ci_is_bioequivalent():
    result = analyze(DATA, parameters=["AUClast"])[0]

    assert result.n_subjects == 8
    assert 100.0 < result.ci.point_estimate < 104.0
    assert result.ci.lower > 80.0
    assert result.ci.upper < 125.0
    assert result.conclusion == "Bioequivalent"


def test_report_file_is_written(tmp_path):
    report = tmp_path / "report.md"

    analyze(DATA, parameters=["AUClast"], report=report)

    assert report.exists()
    assert "BioEqPy Bioequivalence Report" in report.read_text(encoding="utf-8")

