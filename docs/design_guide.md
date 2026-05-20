# Design Guide

## Supported in v0.1

BioEqPy currently supports the standard two-treatment, two-period, two-sequence crossover design:

| Sequence | Period 1 | Period 2 |
|----------|----------|----------|
| TR | Test | Reference |
| RT | Reference | Test |

The dataset must include one row per subject-period observation. Each complete subject must receive both `T` and `R`.

## Required Columns

| Column | Meaning | Example |
|--------|---------|---------|
| `subject` | Subject identifier | `1`, `S001` |
| `sequence` | Assigned treatment sequence | `TR`, `RT` |
| `period` | Period number | `1`, `2` |
| `treatment` | Observed treatment | `T`, `R` |

PK parameter columns must be numeric and positive because the statistical analysis uses natural-log transformed values.

## Analysis Model

The v0.1 2x2 model estimates the treatment effect from subject-level period differences:

```text
period1 - period2 = period_effect + sequence_sign * treatment_difference
```

where:

- `sequence_sign = +1` for `TR`
- `sequence_sign = -1` for `RT`
- `treatment_difference = ln(Test) - ln(Reference)`

The geometric mean ratio is:

```text
GMR = exp(treatment_difference) * 100
```

The 90% CI is:

```text
exp(treatment_difference +/- t(0.95, df) * SE) * 100
```

## Current Limits

The source table includes all five regulatory rows with fully computed sums of squares. The sequence row uses a split-plot F-test with the Subject-within-Sequence mean square as its denominator. SAS decimal-level cross-validation remains outstanding.

## Planned Designs

- 2x2x3 partial replicate.
- 2x2x4 full replicate.
- Williams 2x3x3.
- 4x4 Latin square.
- Parallel group analysis beyond basic detection.

