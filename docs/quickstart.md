# Quickstart

Install the package in editable mode from the project root:

```bash
pip install -e ".[dev]"
```

Run the bundled 2x2 crossover example:

```bash
bioeqpy analyze --data tests/reference_data/example_2x2.csv --parameters AUClast Cmax
```

Estimate a standard 2x2 crossover sample size:

```bash
bioeqpy power --design 2x2 --gmr 0.95 --cv 0.25 --power 0.80
```

The v0.1 implementation is intentionally scoped to the standard 2x2 crossover path. Treat replicate designs, RSABE, and SAS/R decimal-level validation as roadmap items until validated reference datasets are committed.

## Python Example

```python
from bioeqpy import analyze
from bioeqpy.reporting import format_ci_table

results = analyze("tests/reference_data/example_2x2.csv")
print(format_ci_table(results).to_string(index=False))
```

## Documentation Index

- `api_reference.md`: public Python and CLI API.
- `design_guide.md`: supported study design and modeling convention.
- `regulatory_mapping.md`: implemented vs planned regulatory requirements.
- `validation_report.md`: current verification status and validation gaps.
- `roadmap.md`: staged build plan from v0.1 to v1.0.

