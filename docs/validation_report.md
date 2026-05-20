# Validation Report

## Current Verification

The v0.1 test suite passes locally:

```text
18 passed
```

Command used:

```bash
uv run --with pytest --with pandas --with numpy --with scipy --with statsmodels pytest -q
```

## Covered by Tests

- Dataset loading from CSV.
- Numeric PK parameter inference.
- Missing parameter validation.
- 2x2 crossover design detection.
- End-to-end analysis across AUClast, AUCinf, and Cmax.
- 90% CI pass/fail behavior for the bundled example dataset.
- ANOVA source table: no NaN SS values for any parameter.
- ANOVA source table: total DF equals n_observations minus 1.
- ANOVA source table: SS values partition the total SS exactly.
- Sequence F-ratio uses Subject within Sequence as denominator, not Residual.
- Residual MS in source table matches ANOVAResult.residual_ms.
- CI point estimate and limits are unaffected by ANOVA SS computation.
- HTML report writing.
- CLI `analyze`.
- CLI `power`.
- Approximate 2x2 sample-size and power helpers.

## Reference Dataset

The bundled `tests/reference_data/example_2x2.csv` is a synthetic 8-subject 2x2 crossover dataset. It is useful for smoke testing and API demonstration, but it is not a regulatory validation dataset.

## Validation Gaps

Before BioEqPy can be described as regulatorily validated, these files should be added:

- FDA or textbook 2x2 worked example dataset.
- SAS PROC GLM output for the same 2x2 dataset.
- R PowerTOST sample-size reference table.
- Unbalanced 2x2 dropout dataset with expected CI.
- Replicate-design dataset with known within-subject reference variance.

## Next Validation Target

The immediate next target is a decimal-exact 2x2 comparison:

1. Commit a known 2x2 reference dataset from an FDA guidance example or textbook.
2. Commit SAS PROC GLM output as JSON.
3. Assert SS, MS, DF, F, and p-values to the documented precision.
