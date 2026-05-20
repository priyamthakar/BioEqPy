"""Linear-model routines for supported BE designs."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from bioeqpy.core.types import ANOVAResult, BEStudy
from bioeqpy.designs.crossover_2x2 import validate_2x2_completeness


def fit_anova(study: BEStudy) -> ANOVAResult:
    """Fit the appropriate ANOVA model for a study."""
    if study.design.name == "2x2":
        return fit_2x2_anova(study)
    elif study.design.name == "2x2x4":
        return fit_2x2x4_anova(study)
    else:
        raise NotImplementedError(
            f"ANOVA is not yet implemented for design '{study.design.name}'. "
            "Supported designs: 2x2, 2x2x4."
        )


def fit_2x2_anova(study: BEStudy) -> ANOVAResult:
    """Fit a transparent 2x2 crossover ANOVA using subject-level contrasts.

    Treatment estimate and SE come from the period-difference contrast model,
    which is numerically exact and avoids category-ordering footguns.
    Sums of squares for all source table rows are computed from the standard
    between-subject / within-subject decomposition (Type I = Type III for balanced
    2x2 designs).
    """
    frame = validate_2x2_completeness(study)
    pivot = frame.pivot(index="subject", columns="period", values="value")
    subject_sequences = frame.groupby("subject")["sequence"].first().astype(str).str.upper()
    if set(pivot.columns) != {1, 2}:
        raise ValueError("2x2 crossover analysis requires periods 1 and 2.")

    period_difference = pivot[1] - pivot[2]
    sign = subject_sequences.map({"TR": 1.0, "RT": -1.0})
    if sign.isna().any():
        raise ValueError("2x2 crossover analysis requires TR and RT sequence labels.")

    # --- Contrast model for treatment_diff and se_diff ---
    # D_i = beta_0 + sign_i * beta_1 + eps_i
    # beta_1 = mean treatment diff (T minus R), scaled to ln-space
    x = np.column_stack([np.ones(len(period_difference)), sign.to_numpy(dtype=float)])
    y = period_difference.to_numpy(dtype=float)
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    fitted_contrast = pd.Series(x @ beta, index=period_difference.index, name="fitted_contrast")
    contrast_residuals = pd.Series(y - fitted_contrast.to_numpy(), index=period_difference.index)

    residual_df = float(len(y) - x.shape[1])
    sse_contrast = float(np.sum(np.square(contrast_residuals.to_numpy())))
    # Divide by 2: D_i has variance 2*sigma^2; residual_ms estimates sigma^2
    residual_ms = sse_contrast / residual_df / 2.0
    cov_beta = (sse_contrast / residual_df) * np.linalg.inv(x.T @ x)
    treatment_diff = float(beta[1])
    se_diff = float(np.sqrt(cov_beta[1, 1]))

    fitted, residuals = _fit_full_fixed_effects(frame)
    source_table = _build_source_table(frame, residual_ms, residual_df, treatment_diff, se_diff)

    return ANOVAResult(
        source_table=source_table,
        residual_ms=residual_ms,
        residual_df=residual_df,
        treatment_diff=treatment_diff,
        se_diff=se_diff,
        residuals=residuals,
        fitted=fitted,
    )


def fit_2x2x4_anova(study: BEStudy) -> ANOVAResult:
    """Fit ANOVA for 2x2x4 full replicate crossover.

    Estimates treatment difference and within-subject variance components
    for both Test and Reference treatments.
    """
    from bioeqpy.designs.replicate_2x2x4 import validate_2x2x4_completeness
    from bioeqpy.anova.replicate import estimate_within_subject_variances

    frame = validate_2x2x4_completeness(study)

    sWR, sWT, MSW_R, MSW_T = estimate_within_subject_variances(study)

    # Treatment difference: for each subject compute mean_T - mean_R
    subj_diffs = []
    for subject, subj_frame in frame.groupby("subject"):
        t_mean = subj_frame[subj_frame["treatment"] == "T"]["value"].mean()
        r_mean = subj_frame[subj_frame["treatment"] == "R"]["value"].mean()
        subj_diffs.append(t_mean - r_mean)

    treatment_diff = float(np.mean(subj_diffs))

    # SE of treatment difference: sqrt(pooled_MSW / n)
    pooled_MSW = (MSW_R + MSW_T) / 2.0
    n_subjects = len(subj_diffs)
    # For 2x2x4, df_residual = n - 2 (from within-subject contrasts)
    residual_df = float(n_subjects - 2)
    residual_ms = pooled_MSW
    se_diff = float(np.sqrt(pooled_MSW / n_subjects))

    # Treatment F-statistic, SS, and p-value for source table
    treatment_f = (treatment_diff / se_diff) ** 2 if se_diff > 0 else 0.0
    treatment_ss = treatment_f * pooled_MSW
    treatment_p = float(stats.f.sf(treatment_f, 1, residual_df)) if se_diff > 0 else 1.0

    # Build a minimal source table
    source_table = pd.DataFrame([
        {
            "Source": "Treatment",
            "DF": 1,
            "SS": treatment_ss,
            "MS": treatment_ss,
            "F": treatment_f,
            "p-value": treatment_p,
        },
        {
            "Source": "Residual",
            "DF": residual_df,
            "SS": residual_df * residual_ms,
            "MS": residual_ms,
            "F": np.nan,
            "p-value": np.nan,
        },
        {
            "Source": "Within-subject R",
            "DF": n_subjects,
            "SS": n_subjects * MSW_R,
            "MS": MSW_R,
            "F": np.nan,
            "p-value": np.nan,
        },
        {
            "Source": "Within-subject T",
            "DF": n_subjects,
            "SS": n_subjects * MSW_T,
            "MS": MSW_T,
            "F": np.nan,
            "p-value": np.nan,
        },
    ])

    # Residuals from full fixed-effects model
    fitted, residuals = _fit_full_fixed_effects(frame)

    return ANOVAResult(
        source_table=source_table,
        residual_ms=residual_ms,
        residual_df=residual_df,
        treatment_diff=treatment_diff,
        se_diff=se_diff,
        residuals=residuals,
        fitted=fitted,
        swr=sWR,
        swt=sWT,
    )


def _compute_between_subject_ss(frame: pd.DataFrame) -> tuple[float, float, int, int]:
    """Compute Sequence and Subject(Sequence) SS via between-subject decomposition.

    SS_sequence  = Σ_j  n_j × (ȳ_j − ȳ)²   (one-way SS across sequence groups)
    SS_sub_seq   = SS_between − SS_sequence
    SS_between   = Σ_i  n_i × (ȳ_i − ȳ)²   (total between-subject variation)

    Returns (ss_sequence, ss_sub_seq, df_sequence, df_sub_seq).
    For a balanced 2x2 design, Type I equals Type III for these effects.
    """
    grand_mean = float(frame["value"].mean())
    n_subjects = int(frame["subject"].nunique())
    n_sequences = int(frame["sequence"].nunique())

    seq_agg = frame.groupby("sequence")["value"].agg(["mean", "count"])
    ss_seq = float((seq_agg["count"] * (seq_agg["mean"] - grand_mean) ** 2).sum())

    subj_agg = frame.groupby("subject")["value"].agg(["mean", "count"])
    ss_between = float((subj_agg["count"] * (subj_agg["mean"] - grand_mean) ** 2).sum())

    ss_sub_seq = ss_between - ss_seq
    return ss_seq, ss_sub_seq, n_sequences - 1, n_subjects - n_sequences


def _fit_full_fixed_effects(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Fit value ~ subject + period + treatment for residual diagnostics."""
    design = pd.get_dummies(
        frame[["subject", "period", "treatment"]].astype(str),
        columns=["subject", "period", "treatment"],
        drop_first=True,
        dtype=float,
    )
    design.insert(0, "intercept", 1.0)
    x = design.to_numpy(dtype=float)
    y = frame["value"].to_numpy(dtype=float)
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    fitted = pd.Series(x @ beta, index=frame.index, name="fitted")
    residuals = pd.Series(y - fitted.to_numpy(), index=frame.index, name="residual")
    return fitted, residuals


def _effect_pvalue(effect: float, se: float, df: float) -> float:
    """One-df F-test p-value: F = (effect/se)^2, p = Pr(F_{1,df} > F)."""
    if se == 0:
        return 1.0
    f_value = (effect / se) ** 2
    return float(stats.f.sf(f_value, 1, df))


def _build_source_table(
    frame: pd.DataFrame,
    residual_ms: float,
    residual_df: float,
    treatment_diff: float,
    se_diff: float,
) -> pd.DataFrame:
    """Build the regulatory ANOVA source table with all SS values populated.

    Split-plot F-test structure (per FDA 2001 / EMA 2010):
      Sequence         → F = MS_seq / MS_sub_seq
      Subject(Seq)     → F = MS_sub_seq / MS_residual
      Period, Treatment → F = MS_effect / MS_residual
    """
    n_subjects = int(frame["subject"].nunique())

    # --- Between-subject effects ---
    ss_seq, ss_sub_seq, df_seq, df_sub_seq = _compute_between_subject_ss(frame)
    ms_seq = ss_seq / df_seq if df_seq > 0 else np.nan
    ms_sub_seq = ss_sub_seq / df_sub_seq if df_sub_seq > 0 else np.nan
    # Sequence: split-plot denominator is Subject(Sequence)
    if not np.isnan(ms_sub_seq) and ms_sub_seq > 0:
        f_seq = float(ms_seq / ms_sub_seq)
        p_seq = float(stats.f.sf(f_seq, df_seq, df_sub_seq))
    else:
        f_seq = np.nan
        p_seq = np.nan
    # Subject(Sequence): tested against residual
    if not np.isnan(ms_sub_seq) and residual_ms > 0:
        f_sub_seq = float(ms_sub_seq / residual_ms)
        p_sub_seq = float(stats.f.sf(f_sub_seq, df_sub_seq, residual_df))
    else:
        f_sub_seq = np.nan
        p_sub_seq = np.nan

    # --- Within-subject effects ---
    period_means = frame.groupby("period")["value"].mean()
    period_effect = float(period_means.get(1, 0.0) - period_means.get(2, 0.0))
    period_se = float(np.sqrt(max(residual_ms, 0.0) * 2.0 / n_subjects))
    period_f = (period_effect / period_se) ** 2 if period_se > 0 else 0.0
    period_ss = period_f * residual_ms

    treatment_f = (treatment_diff / se_diff) ** 2 if se_diff > 0 else 0.0
    treatment_ss = treatment_f * residual_ms
    residual_ss = residual_ms * residual_df

    rows = [
        {
            "Source": "Sequence",
            "DF": df_seq,
            "SS": ss_seq,
            "MS": ms_seq,
            "F": f_seq,
            "p-value": p_seq,
        },
        {
            "Source": "Subject within Sequence",
            "DF": df_sub_seq,
            "SS": ss_sub_seq,
            "MS": ms_sub_seq,
            "F": f_sub_seq,
            "p-value": p_sub_seq,
        },
        {
            "Source": "Period",
            "DF": 1,
            "SS": period_ss,
            "MS": period_ss,
            "F": period_f,
            "p-value": _effect_pvalue(period_effect, period_se, residual_df),
        },
        {
            "Source": "Treatment",
            "DF": 1,
            "SS": treatment_ss,
            "MS": treatment_ss,
            "F": treatment_f,
            "p-value": _effect_pvalue(treatment_diff, se_diff, residual_df),
        },
        {
            "Source": "Residual",
            "DF": residual_df,
            "SS": residual_ss,
            "MS": residual_ms,
            "F": np.nan,
            "p-value": np.nan,
        },
    ]
    return pd.DataFrame(rows)
