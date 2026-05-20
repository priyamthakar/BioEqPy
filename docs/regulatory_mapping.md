# Regulatory Mapping

BioEqPy is structured around standard FDA and EMA average bioequivalence requirements. This document states what the v0.1 code implements and what remains planned.

## Implemented

| Requirement | v0.1 Status |
|-------------|-------------|
| Natural-log transformation of AUC and Cmax | Implemented |
| 90% CI for geometric mean ratio | Implemented |
| Standard 80.00-125.00% acceptance limits | Implemented |
| User-specified acceptance limits | Implemented through `analyze(..., acceptance=(lower, upper))` |
| Standard 2x2 crossover treatment comparison | Implemented |
| Basic residual normality diagnostic | Implemented |
| Basic outlier flag using standardized residuals | Implemented |
| Sample-size planning for 2x2 crossover | Implemented as approximate closed-form calculation |

## Partially Implemented

| Requirement | Current Status |
|-------------|----------------|
| ANOVA source table | All rows fully computed; SAS decimal-level cross-validation still outstanding |
| Period effect test | Included as an approximate fixed-effect test |
| Sequence effect test | Implemented with split-plot F-test |
| Subject-within-Sequence row | Implemented; SS computed from between-subject decomposition |

## Not Yet Implemented

| Requirement | Planned Module |
|-------------|----------------|
| Reference-scaled average bioequivalence | `bioeqpy.tost.scaled` |
| Hyslop upper confidence bound | `bioeqpy.tost.scaled` |
| Replicate within-subject variance components | `bioeqpy.anova` and replicate design modules |
| NTI-specific workflow | `bioeqpy.tost.nti` |
| EMA replicate-specific reporting | `bioeqpy.reporting` |
| SAS PROC GLM decimal validation | `tests/test_anova_vs_sas.py` |
| PowerTOST exact sample-size validation | `tests/test_power.py` with committed reference values |

## Practical Interpretation

Use the current package as a working Phase 1 implementation for standard 2x2 crossover examples, portfolio demonstration, and API iteration. Do not claim regulatory production equivalence until the SAS/R validation datasets and decimal comparison tests are committed.

