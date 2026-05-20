# BioEqPy Project Context

## Current State

Phase 1 is feature-complete and all 18 tests pass (1.53s). The package is installable,
CLI-ready, and generates a professional HTML report.

## Implemented

- Package metadata in `pyproject.toml` with optional `[reports]` extra for jinja2.
- Public API: `bioeqpy.analyze`.
- Core dataclasses: `BEStudy`, `DesignSpec`, `ANOVAResult`, `CIResult`, `BEResult`.
- CSV/Excel dataset loading and validation.
- 2x2 crossover design detection.
- Subject-level contrast model for `treatment_diff` and `se_diff` (numerically exact).
- **Complete ANOVA source table** (no NaN placeholders):
  - Sequence SS/MS/F via between-subject one-way decomposition.
  - Subject(Sequence) SS/MS/F as residual of between-subject variation.
  - Split-plot F-test structure: Sequence tested against Subject(Sequence) MS.
  - Period and Treatment tested against Residual MS.
  - SS partition verified: Σ SS = SS_total (test_anova_table.py).
- Standard TOST 90% CI and pass/fail decision.
- Summary statistics and compact BE assessment table.
- **Aesthetic HTML report** (`analyze(..., report="report.html")`):
  - Self-contained single file (inline CSS + SVG, no external JS).
  - Pharma blue (`#1A4B8C`) design system with Inter + JetBrains Mono typography.
  - Per-parameter forest plots (inline SVG, acceptance zone shaded, diamond marker).
  - Split ANOVA source table with Treatment row highlighted.
  - Diagnostics grid (residual MS/df, Shapiro-Wilk p, outlier count).
  - Print-optimised with `@media print` CSS.
- Approximate 2x2 sample-size and power helpers.
- CLI commands: `bioeqpy analyze` and `bioeqpy power`.
- Example reference dataset and 18 tests.
- Documentation set under `docs/`.
- `example_report.html` in project root (generated from example data).

## Test Suite

```
18 passed in 1.53s
```

Key structural tests added in `tests/test_anova_table.py`:
- SS has no NaN values.
- DF sums to n_obs − 1 = 15.
- SS_table == SS_total (partition correctness to 1e-8).
- Sequence F = MS_seq / MS_sub_seq (split-plot denominator).
- Residual MS in table matches `ANOVAResult.residual_ms`.

## Not Yet Implemented

- Replicate designs (2x2x3 partial, 2x2x4 full, Williams 2x3x3, Latin 4x4).
- RSABE / Hyslop scaled criterion for HVDs (sWR > 0.294).
- NTI drug tightened limits (90.00–111.11%).
- Carryover effect testing.
- HTML report: power curve plot, Q-Q plot SVG.
- HTML report: PDF export via WeasyPrint.
- OpenPKFlow NCA bridge.
- SAS/R decimal-exact cross-validation (documented as outstanding).
- CI/CD and PyPI release workflow.

## Next Technical Steps (priority order)

1. **External validation**: Run the AUClast ANOVA against R `lm()` or SAS PROC GLM;
   record expected SS values as JSON sidecar for `tests/test_anova_vs_reference.py`.
2. **Replicate designs**: Implement 2x2x4 full replicate with REML within-subject variance
   components; enables RSABE for HVDs.
3. **RSABE criterion**: Hyslop linearised criterion, FDA 2022 Progesterone methodology.
4. **Improved diagnostics**: Cook's distance, leverage, Q-Q data for the HTML report.
5. **Report extras**: power curve SVG panel, Q-Q residual plot, WeasyPrint PDF export.
