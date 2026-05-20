# Roadmap

## v0.1 Completed

- Package scaffold and core dataclasses.
- CSV and Excel loading with column inference.
- 2x2 crossover design detection and validation.
- Log-transformation of raw PK parameters.
- Full split-plot ANOVA source table: Sequence, Subject within Sequence, Period, Treatment, Residual — all SS values populated, no placeholders.
- Sequence tested against Subject within Sequence as denominator per FDA 2001 and EMA 2010.
- Standard TOST 90% CI on geometric mean ratio.
- 80-125% acceptance limits with pass/fail conclusion.
- HTML report with inline SVG forest plots via Jinja2 renderer.
- Basic diagnostics: normality check and outlier flagging.
- Sample-size and power estimation for 2x2 crossover.
- CLI: `bioeqpy analyze` and `bioeqpy power`.
- Synthetic example dataset for smoke testing.
- 18 passing tests, including structural ANOVA table validation.
- User-facing documentation suite.

## v0.2 Target

- Partial and full replicate design loaders.
- Narrow therapeutic index handling with tightened acceptance limits.
- Carryover test for 2x2 crossover.
- Within-subject variance estimates for Reference and Test.

## v0.3 Target

- RSABE for highly variable drugs, with Hyslop criterion.
- Full replicate ANOVA source table.
- Expanded CLI options for alpha and acceptance limit overrides.

## v0.4 Target

- PDF report output.
- Power curve SVG embedded in HTML report.
- OpenPKFlow bridge for PK import.

## v1.0 Target

- SAS and R decimal validation for every advertised method.
- 200+ tests.
- PyPI release.
- JOSS or SoftwareX-ready validation narrative.
