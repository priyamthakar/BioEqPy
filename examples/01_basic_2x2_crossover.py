"""Basic 2x2 crossover bioequivalence analysis."""

from pathlib import Path

from bioeqpy import analyze
from bioeqpy.reporting import format_ci_table


data_path = Path(__file__).resolve().parents[1] / "tests" / "reference_data" / "example_2x2.csv"
results = analyze(data_path, parameters=["AUClast", "Cmax"])
print(format_ci_table(results).to_string(index=False))

