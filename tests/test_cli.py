from pathlib import Path

from bioeqpy.cli.main import main


DATA = Path(__file__).parent / "reference_data" / "example_2x2.csv"


def test_cli_analyze(capsys):
    status = main(["analyze", "--data", str(DATA), "--parameters", "AUClast"])
    captured = capsys.readouterr()

    assert status == 0
    assert "AUClast" in captured.out
    assert "90% CI" in captured.out


def test_cli_power(capsys):
    status = main(["power", "--cv", "0.25", "--gmr", "0.95", "--power", "0.8"])
    captured = capsys.readouterr()

    assert status == 0
    assert "Required N:" in captured.out

