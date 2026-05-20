# API Reference

This document covers the v0.1 public API. Functions not listed here should be treated as internal implementation details until the package reaches a broader validation milestone.

## `bioeqpy.analyze`

```python
from bioeqpy import analyze

results = analyze(
    data="tests/reference_data/example_2x2.csv",
    parameters=["AUClast", "Cmax"],
    alpha=0.05,
    acceptance=(80.0, 125.0),
    method="standard",
    report="be_report.md",
)
```

Runs the standard 2x2 crossover workflow for one or more PK parameters.

Inputs:

- `data`: CSV path, Excel path, or a `pandas.DataFrame`.
- `parameters`: optional list of PK parameter columns. If omitted, all numeric non-design columns are analyzed.
- `alpha`: one-sided alpha for TOST. The regulatory default is `0.05`, producing a 90% CI.
- `acceptance`: percentage-scale acceptance limits. The default is `(80.0, 125.0)`.
- `method`: method label stored in the result. In v0.1, only standard average bioequivalence is implemented.
- `report`: optional report path. Use a `.html` extension for a styled HTML report; otherwise a plain Markdown file is written.

Returns a list of `BEResult` objects, one per parameter.

## Data Containers

`BEStudy` stores one log-transformed PK parameter dataset.

`BEResult` stores the complete analysis result for one parameter, including summary statistics, ANOVA output, CI output, diagnostics, conclusion, and regulatory basis.

`ANOVAResult` stores the treatment estimate, standard error, residual degrees of freedom, residual mean square, source table, residuals, and fitted values.

`CIResult` stores the geometric mean ratio, CI limits, acceptance limits, pass/fail flag, and method label.

## Loading Data

```python
from bioeqpy.io import load_dataset, load_study

data = load_dataset("be_data.csv")
study = load_study(data, parameter="AUClast")
```

`load_dataset` accepts `.csv`, `.xlsx`, and `.xls`. `load_study` validates the selected PK parameter, applies natural-log transformation, and attaches the detected design.

## Power Helpers

```python
from bioeqpy.power import achieved_power, sample_size, power_curve

n = sample_size(cv=0.25, gmr=0.95, power=0.80)
power = achieved_power(n=24, cv=0.25, gmr=0.95)
curve = power_curve(range(12, 60, 2), cv=0.25)
```

The v0.1 power implementation is an approximate closed-form TOST calculation for balanced 2x2 crossover studies.

## CLI

```bash
bioeqpy analyze --data tests/reference_data/example_2x2.csv --parameters AUClast Cmax
bioeqpy analyze --data tests/reference_data/example_2x2.csv --report be_report.md
bioeqpy power --design 2x2 --gmr 0.95 --cv 0.25 --power 0.80
bioeqpy --version
```

