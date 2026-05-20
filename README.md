# BioEqPy

Python bioequivalence statistics library. No SAS or R required.

BioEqPy implements the standard 2x2 crossover BE workflow: log-transform raw PK parameters, run a split-plot ANOVA, compute the 90% confidence interval via TOST, and emit a professional HTML report with inline SVG forest plots. The library is pip-installable and exposes both a Python API and a CLI.

---

## Install

```bash
pip install -e .
```

With HTML reports:

```bash
pip install -e ".[reports]"
```

Runtime dependencies: numpy, pandas, scipy, statsmodels. Jinja2 is needed only for HTML output.

---

## Quickstart

```python
from bioeqpy import analyze

# Returns a list of BEResult objects, one per PK parameter
results = analyze("study.csv")

for r in results:
    print(r.parameter_name, r.ci.point_estimate, r.ci.lower, r.ci.upper)

# HTML report with inline SVG forest plots
analyze("study.csv", report="be_report.html")
```

### CSV format

Required columns: `subject`, `sequence`, `period`, `treatment`. Any additional numeric columns are treated as PK parameters.

```csv
subject,sequence,period,treatment,AUClast,AUCinf,Cmax
1,TR,1,T,520.3,580.1,45.2
1,TR,2,R,498.7,555.0,42.8
2,RT,1,R,601.2,670.5,50.1
2,RT,2,T,585.0,648.3,48.7
```

Valid sequence labels: `TR` and `RT`. Periods: `1` and `2`. Treatments: `T` and `R`.

---

## CLI

```bash
# Run analysis and print the CI table
bioeqpy analyze --data study.csv --parameters AUClast Cmax

# Write an HTML report
bioeqpy analyze --data study.csv --report be_report.html

# Estimate required sample size
bioeqpy power --cv 0.25 --gmr 0.95 --power 0.8
```

---

## Features

| Feature | Status |
|---|---|
| Standard 2x2 crossover analysis | ✓ |
| ANOVA source table with split-plot F-tests | ✓ |
| 90% CI on geometric mean ratio via TOST | ✓ |
| 80-125% acceptance limits | ✓ |
| HTML report with inline SVG forest plots | ✓ |
| CLI: `bioeqpy analyze` and `bioeqpy power` | ✓ |
| Sample-size and power estimation for 2x2 crossover | ✓ |
| Replicate designs | planned |
| Reference-scaled ABE for highly variable drugs | planned |
| Narrow therapeutic index handling | planned |

---

## Design

- All PK parameters are log-transformed before analysis. The geometric mean ratio and its 90% CI are back-transformed and expressed as percentages.
- The ANOVA follows the split-plot structure specified in FDA 2001 and EMA 2010 guidance: Sequence is tested against Subject within Sequence; Period and Treatment are tested against the Residual.
- The TOST procedure checks whether both one-sided t-tests reject at the 0.05 level, which is algebraically equivalent to the 90% CI falling within 80-125%.
- The forest plot SVG is self-contained and embedded inline so the HTML report has no external dependencies.

---

## Status

Phase 1 complete. 18 tests passing. Not yet validated against SAS or R reference output to decimal precision.

```
18 passed
```

---

## Docs

- [Quickstart](docs/quickstart.md)
- [API Reference](docs/api_reference.md)
- [Design Guide](docs/design_guide.md)
- [Regulatory Mapping](docs/regulatory_mapping.md)
- [Validation Report](docs/validation_report.md)
- [Roadmap](docs/roadmap.md)

---

## License

MIT. Author: Priyam T.
