"""High-level BioEqPy analysis orchestration."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from bioeqpy.anova import fit_anova
from bioeqpy.core.constants import FDA_2001_STANDARD_BASIS, FDA_2011_RSABE_BASIS, EMA_2018_ABEL_BASIS, NTI_ACCEPTANCE, STANDARD_ACCEPTANCE
from bioeqpy.core.types import BEResult
from bioeqpy.diagnostics import check_assumptions
from bioeqpy.io.loaders import infer_parameters, load_dataset, load_study
from bioeqpy.reporting.tables import format_ci_table, summary_statistics
from bioeqpy.tost import compute_ci


def analyze(
    data: str | Path | pd.DataFrame,
    parameters: list[str] | None = None,
    alpha: float = 0.05,
    acceptance: tuple[float, float] = STANDARD_ACCEPTANCE,
    method: str = "standard",
    report: str | Path | None = None,
) -> list[BEResult]:
    """Run BE analysis for one or more PK parameters."""
    dataset = load_dataset(data) if not isinstance(data, pd.DataFrame) else data.copy()
    selected = parameters or infer_parameters(dataset)
    if not selected:
        raise ValueError("No numeric PK parameter columns were found.")

    results: list[BEResult] = []
    for parameter in selected:
        study = load_study(dataset, parameter=parameter)
        anova_result = fit_anova(study)
        if method == "nti":
            from bioeqpy.tost.nti import compute_nti_ci
            ci_result = compute_nti_ci(anova_result, alpha=alpha)
        elif method == "abel":
            from bioeqpy.tost.rsabe import compute_abel_ci
            ci_result = compute_abel_ci(anova_result, alpha=alpha)
        elif method == "rsabe":
            from bioeqpy.tost.rsabe import compute_rsabe_ci
            ci_result = compute_rsabe_ci(anova_result, alpha=alpha)
        else:
            ci_result = compute_ci(anova_result, alpha=alpha, acceptance=acceptance, method=method)
        diagnostics = check_assumptions(study, anova_result)
        basis_map = {
            "nti": FDA_2001_STANDARD_BASIS,
            "abel": EMA_2018_ABEL_BASIS,
            "rsabe": FDA_2011_RSABE_BASIS,
        }
        regulatory_basis = basis_map.get(method, FDA_2001_STANDARD_BASIS)
        results.append(
            BEResult(
                parameter_name=parameter,
                n_subjects=int(pd.Series(study.subjects).nunique()),
                design=study.design,
                summary_stats=summary_statistics(study),
                anova=anova_result,
                ci=ci_result,
                diagnostics=diagnostics,
                conclusion="Bioequivalent" if ci_result.passed else "Not bioequivalent",
                regulatory_basis=regulatory_basis,
            )
        )

    if report is not None:
        _write_text_report(results, Path(report))
    return results


def _write_text_report(results: list[BEResult], path: Path) -> None:
    """Write an HTML report (if .html) or a plain-text Markdown fallback."""
    if path.suffix.lower() in {".html", ".htm"}:
        try:
            from bioeqpy.reporting.renderer import render_html_report
            render_html_report(results, path)
            return
        except ImportError:
            pass  # fall through to text if jinja2 unavailable

    lines = ["# BioEqPy Bioequivalence Report", ""]
    lines.append("```")
    lines.append(format_ci_table(results).to_string(index=False))
    lines.append("```")
    lines.append("")
    for result in results:
        lines.extend(
            [
                f"## {result.parameter_name}",
                "",
                "### Summary Statistics",
                "```",
                result.summary_stats.to_string(index=False),
                "```",
                "",
                "### ANOVA",
                "```",
                result.anova.source_table.to_string(index=False),
                "```",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
