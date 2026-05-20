"""Jinja2-based HTML report renderer for bioequivalence analyses."""

from __future__ import annotations

import datetime
import math
from pathlib import Path
from typing import TYPE_CHECKING

try:
    import jinja2
    _JINJA2_AVAILABLE = True
except ImportError:
    _JINJA2_AVAILABLE = False

if TYPE_CHECKING:
    from bioeqpy.core.types import BEResult

_VERSION = "0.1.0"


def render_html_report(
    results: list[BEResult],
    path: Path,
    title: str = "Bioequivalence Statistical Report",
) -> None:
    """Render a self-contained HTML bioequivalence report."""
    if not _JINJA2_AVAILABLE:
        raise ImportError(
            "Jinja2 is required for HTML reports. "
            "Install it: pip install 'bioeqpy[reports]' or pip install jinja2"
        )

    template_dir = Path(__file__).parent / "templates"
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_dir)),
        autoescape=jinja2.select_autoescape(["html"]),
        undefined=jinja2.StrictUndefined,
    )
    env.filters["fmt_float"] = _fmt_float
    env.filters["fmt_pval"] = _fmt_pval
    env.filters["fmt_ss"] = _fmt_ss

    template = env.get_template("be_report.html")

    prepared = [_prepare_result(r) for r in results]
    all_passed = all(r["passed"] for r in prepared)

    html = template.render(
        title=title,
        results=prepared,
        generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC"),
        version=_VERSION,
        all_passed=all_passed,
        n_parameters=len(results),
    )
    path.write_text(html, encoding="utf-8")


# ─── Jinja2 filter helpers ───────────────────────────────────────────────────

def _fmt_float(value: object, decimals: int = 4) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "—"
    try:
        return f"{float(value):.{decimals}f}"
    except (TypeError, ValueError):
        return str(value)


def _fmt_pval(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "—"
    try:
        p = float(value)
        if p < 0.001:
            return "<0.001"
        if p < 0.01:
            return f"{p:.4f}"
        return f"{p:.4f}"
    except (TypeError, ValueError):
        return str(value)


def _fmt_ss(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "—"
    try:
        v = float(value)
        if abs(v) < 0.001:
            return f"{v:.6f}"
        return f"{v:.5f}"
    except (TypeError, ValueError):
        return str(value)


# ─── Data preparation ─────────────────────────────────────────────────────────

def _prepare_result(result: BEResult) -> dict:
    anova_rows = result.anova.source_table.to_dict("records")
    summary_rows = result.summary_stats.to_dict("records")

    forest_svg = _build_forest_svg(
        point=result.ci.point_estimate,
        lower=result.ci.lower,
        upper=result.ci.upper,
        lo_limit=result.ci.acceptance_lower,
        hi_limit=result.ci.acceptance_upper,
        passed=result.ci.passed,
    )

    power_svg = _build_power_curve_svg(
        n_subjects=result.n_subjects,
        residual_ms=result.anova.residual_ms,
        residual_df=result.anova.residual_df,
        treatment_diff=result.anova.treatment_diff,
    )

    qq_svg = _build_qq_svg(
        residuals_list=result.anova.residuals.dropna().tolist(),
    )

    warnings = result.diagnostics.get("warnings", [])
    normality_p = result.diagnostics.get("normality_p_value")
    outliers = result.diagnostics.get("outlier_indices", [])

    return {
        "parameter_name": result.parameter_name,
        "n_subjects": result.n_subjects,
        "design_name": result.design.name,
        "sequences": result.design.sequences,
        "point_estimate": round(result.ci.point_estimate, 2),
        "ci_lower": round(result.ci.lower, 2),
        "ci_upper": round(result.ci.upper, 2),
        "acceptance_lower": result.ci.acceptance_lower,
        "acceptance_upper": result.ci.acceptance_upper,
        "passed": result.ci.passed,
        "conclusion": result.conclusion,
        "regulatory_basis": result.regulatory_basis,
        "residual_ms": result.anova.residual_ms,
        "residual_df": int(result.anova.residual_df),
        "anova_rows": anova_rows,
        "summary_rows": summary_rows,
        "normality_p": normality_p,
        "outliers": outliers,
        "warnings": warnings,
        "forest_svg": forest_svg,
        "power_svg": power_svg,
        "qq_svg": qq_svg,
    }


# ─── Forest plot SVG ─────────────────────────────────────────────────────────

def _build_forest_svg(
    point: float,
    lower: float,
    upper: float,
    lo_limit: float = 80.0,
    hi_limit: float = 125.0,
    passed: bool = True,
    pct_min: float = 68.0,
    pct_max: float = 142.0,
) -> str:
    """Build a self-contained inline SVG forest plot for one CI row."""
    W, H = 560, 88
    PAD_L, PAD_R, PAD_TOP = 8, 8, 6
    AXIS_Y = 54
    BAR_H = 10  # half-height of end-caps
    DIAMOND = 7  # half-size of diamond

    plot_x0 = PAD_L
    plot_x1 = W - PAD_R

    def px(pct: float) -> float:
        return plot_x0 + (pct - pct_min) / (pct_max - pct_min) * (plot_x1 - plot_x0)

    ci_color = "#16A34A" if passed else "#DC2626"
    zone_fill = "#EFF6FF"
    tick_color = "#CBD5E1"
    text_color = "#475569"
    ref_color = "#94A3B8"
    label_color = ci_color

    ticks = [70, 80, 90, 100, 110, 120, 125, 130, 140]

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="{W}" height="{H}" '
        f'style="font-family:\'JetBrains Mono\',\'Courier New\',monospace;font-size:9.5px;'
        f'display:block;max-width:100%;">'
    )

    # Acceptance zone fill
    z0, z1 = px(lo_limit), px(hi_limit)
    parts.append(
        f'<rect x="{z0:.1f}" y="{PAD_TOP}" width="{z1-z0:.1f}" '
        f'height="{AXIS_Y - PAD_TOP + 4}" fill="{zone_fill}" rx="2"/>'
    )

    # Acceptance limit dashed verticals
    for lim in [lo_limit, hi_limit]:
        lx = px(lim)
        parts.append(
            f'<line x1="{lx:.1f}" y1="{PAD_TOP}" x2="{lx:.1f}" '
            f'y2="{AXIS_Y + 4}" stroke="#94A3B8" stroke-width="1" stroke-dasharray="4,3"/>'
        )

    # Reference line at 100%
    r100 = px(100.0)
    parts.append(
        f'<line x1="{r100:.1f}" y1="{PAD_TOP}" x2="{r100:.1f}" '
        f'y2="{AXIS_Y + 4}" stroke="{ref_color}" stroke-width="0.8" stroke-dasharray="5,4" opacity="0.5"/>'
    )

    # Axis baseline
    parts.append(
        f'<line x1="{plot_x0}" y1="{AXIS_Y}" x2="{plot_x1}" y2="{AXIS_Y}" '
        f'stroke="{tick_color}" stroke-width="0.75"/>'
    )

    # Tick marks + labels
    for tick in ticks:
        tx = px(tick)
        is_key = tick in {80, 100, 125}
        t_color = "#1A4B8C" if is_key else text_color
        weight = "600" if is_key else "400"
        parts.append(
            f'<line x1="{tx:.1f}" y1="{AXIS_Y}" x2="{tx:.1f}" '
            f'y2="{AXIS_Y + 4}" stroke="{tick_color}" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{tx:.1f}" y="{AXIS_Y + 15}" text-anchor="middle" '
            f'fill="{t_color}" font-weight="{weight}">{tick}</text>'
        )

    # CI line (thick horizontal)
    lx_ci, ux_ci = px(lower), px(upper)
    parts.append(
        f'<line x1="{lx_ci:.1f}" y1="{AXIS_Y}" x2="{ux_ci:.1f}" y2="{AXIS_Y}" '
        f'stroke="{ci_color}" stroke-width="3.5" stroke-linecap="round" opacity="0.9"/>'
    )

    # End-cap verticals
    for ex in [lx_ci, ux_ci]:
        parts.append(
            f'<line x1="{ex:.1f}" y1="{AXIS_Y - BAR_H}" x2="{ex:.1f}" '
            f'y2="{AXIS_Y + BAR_H}" stroke="{ci_color}" stroke-width="2"/>'
        )

    # Diamond at point estimate
    pdx = px(point)
    d = DIAMOND
    parts.append(
        f'<polygon points="{pdx:.1f},{AXIS_Y - d} {pdx + d:.1f},{AXIS_Y} '
        f'{pdx:.1f},{AXIS_Y + d} {pdx - d:.1f},{AXIS_Y}" fill="{ci_color}"/>'
    )

    # CI label above the bar
    mid_x = (lx_ci + ux_ci) / 2
    label = f"{point:.2f}%  [{lower:.2f} – {upper:.2f}]"
    parts.append(
        f'<text x="{mid_x:.1f}" y="{PAD_TOP + 11}" text-anchor="middle" '
        f'fill="{label_color}" font-weight="700" font-size="11px">{label}</text>'
    )

    # Limit labels
    parts.append(
        f'<text x="{z0:.1f}" y="{H - 2}" text-anchor="middle" '
        f'fill="#94A3B8" font-size="8.5px">{lo_limit:.0f}%</text>'
    )
    parts.append(
        f'<text x="{z1:.1f}" y="{H - 2}" text-anchor="middle" '
        f'fill="#94A3B8" font-size="8.5px">{hi_limit:.0f}%</text>'
    )

    parts.append("</svg>")
    return "\n".join(parts)


# ─── Power curve SVG ──────────────────────────────────────────────────────────

def _build_power_curve_svg(
    n_subjects: int,
    residual_ms: float,
    residual_df: float,
    treatment_diff: float,
    width: int = 320,
    height: int = 140,
) -> str:
    """Build a power curve SVG showing power vs N for the estimated CV."""
    import math
    from scipy import stats as scipy_stats

    # Estimate CV from residual MS: CV ≈ sqrt(exp(MSW) - 1) ≈ sqrt(MSW) for small MSW
    cv = float(math.sqrt(residual_ms))
    gmr = float(math.exp(treatment_diff))

    # Compute power for a range of N values
    n_values = list(range(4, max(n_subjects * 2 + 4, 40), 2))

    def power_at_n(n, cv, gmr, alpha=0.05, lo=0.8, hi=1.25):
        # Approximate TOST power for 2x2 crossover
        df = n - 2
        if df <= 0:
            return 0.0
        se = math.sqrt(2 * cv**2 / n)
        tcrit = scipy_stats.t.ppf(1 - alpha, df)
        delta_lo = (math.log(lo) - math.log(gmr)) / se
        delta_hi = (math.log(hi) - math.log(gmr)) / se
        power = (
            scipy_stats.t.cdf(delta_hi - tcrit, df)
            - scipy_stats.t.cdf(delta_lo + tcrit, df)
        )
        return max(0.0, float(power))

    powers = [power_at_n(n, cv, gmr) for n in n_values]

    # SVG parameters
    W, H = width, height
    PAD_L, PAD_R, PAD_TOP, PAD_BOT = 36, 12, 10, 28
    plot_w = W - PAD_L - PAD_R
    plot_h = H - PAD_TOP - PAD_BOT

    def sx(n_val):
        n_min, n_max = n_values[0], n_values[-1]
        return PAD_L + (n_val - n_min) / (n_max - n_min) * plot_w

    def sy(p_val):
        return PAD_TOP + (1.0 - p_val) * plot_h

    # Build path
    pts = " ".join(f"{sx(n):.1f},{sy(p):.1f}" for n, p in zip(n_values, powers))

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
        f'style="font-family:\'JetBrains Mono\',monospace;font-size:9px;display:block;">',
    ]

    # Background
    lines.append(f'<rect x="{PAD_L}" y="{PAD_TOP}" width="{plot_w}" height="{plot_h}" fill="#F8FAFC" rx="3" stroke="#E2E8F0" stroke-width="0.5"/>')

    # 80% power reference line
    p80_y = sy(0.80)
    lines.append(f'<line x1="{PAD_L}" y1="{p80_y:.1f}" x2="{PAD_L + plot_w}" y2="{p80_y:.1f}" stroke="#16A34A" stroke-width="0.75" stroke-dasharray="4,3" opacity="0.7"/>')
    lines.append(f'<text x="{PAD_L + 3}" y="{p80_y - 3:.1f}" fill="#16A34A" font-size="8px">80%</text>')

    # Current N vertical line
    cur_x = sx(n_subjects)
    lines.append(f'<line x1="{cur_x:.1f}" y1="{PAD_TOP}" x2="{cur_x:.1f}" y2="{PAD_TOP + plot_h}" stroke="#1A4B8C" stroke-width="1" stroke-dasharray="3,2" opacity="0.6"/>')

    # Power curve
    lines.append(f'<polyline points="{pts}" fill="none" stroke="#1A4B8C" stroke-width="2" stroke-linejoin="round"/>')

    # Dot at current N
    cur_power = power_at_n(n_subjects, cv, gmr)
    cur_y = sy(cur_power)
    lines.append(f'<circle cx="{cur_x:.1f}" cy="{cur_y:.1f}" r="4" fill="#1A4B8C"/>')
    lines.append(f'<text x="{cur_x + 5:.1f}" y="{cur_y - 4:.1f}" fill="#1A4B8C" font-weight="600" font-size="9px">n={n_subjects}, {cur_power*100:.0f}%</text>')

    # Y axis ticks
    for p_tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        ty = sy(p_tick)
        lines.append(f'<line x1="{PAD_L - 3}" y1="{ty:.1f}" x2="{PAD_L}" y2="{ty:.1f}" stroke="#CBD5E1" stroke-width="1"/>')
        lines.append(f'<text x="{PAD_L - 5}" y="{ty + 3:.1f}" text-anchor="end" fill="#94A3B8">{int(p_tick*100)}%</text>')

    # X axis ticks (a few)
    n_range = n_values[-1] - n_values[0]
    step = 8 if n_range <= 32 else 16
    for n_tick in range(n_values[0], n_values[-1] + 1, step):
        tx = sx(n_tick)
        lines.append(f'<line x1="{tx:.1f}" y1="{PAD_TOP + plot_h}" x2="{tx:.1f}" y2="{PAD_TOP + plot_h + 3}" stroke="#CBD5E1" stroke-width="1"/>')
        lines.append(f'<text x="{tx:.1f}" y="{PAD_TOP + plot_h + 13}" text-anchor="middle" fill="#94A3B8">{n_tick}</text>')

    # Axis labels
    lines.append(f'<text x="{PAD_L + plot_w // 2}" y="{H - 2}" text-anchor="middle" fill="#64748B" font-size="9px">Sample Size N</text>')
    lines.append(f'<text x="9" y="{PAD_TOP + plot_h // 2}" text-anchor="middle" fill="#64748B" font-size="9px" transform="rotate(-90,9,{PAD_TOP + plot_h // 2})">Power</text>')

    # CV annotation
    lines.append(f'<text x="{PAD_L + plot_w - 2}" y="{PAD_TOP + 11}" text-anchor="end" fill="#94A3B8" font-size="8px">CV ≈ {cv*100:.1f}%</text>')

    lines.append('</svg>')
    return '\n'.join(lines)


# ─── Q-Q residual plot SVG ────────────────────────────────────────────────────

def _build_qq_svg(
    residuals_list: list[float],
    width: int = 280,
    height: int = 180,
) -> str:
    """Build a Q-Q plot SVG for residuals vs theoretical normal quantiles."""
    import math

    if len(residuals_list) < 3:
        return '<svg width="280" height="80"><text x="10" y="40" fill="#94A3B8" font-size="12px">Insufficient residuals for Q-Q plot</text></svg>'

    n = len(residuals_list)
    sorted_res = sorted(residuals_list)
    mean_r = sum(sorted_res) / n
    std_r = math.sqrt(sum((r - mean_r)**2 for r in sorted_res) / (n - 1))

    # Theoretical quantiles using Blom's formula: (i - 0.375) / (n + 0.25)
    def norm_ppf(p):
        # Beasley-Springer-Moro approximation for normal quantile
        if p <= 0 or p >= 1:
            return 0.0
        a = [2.515517, 0.802853, 0.010328]
        b = [1.432788, 0.189269, 0.001308]
        if p > 0.5:
            t = math.sqrt(-2 * math.log(1 - p))
            sign = 1
        else:
            t = math.sqrt(-2 * math.log(p))
            sign = -1
        num = a[0] + a[1]*t + a[2]*t*t
        den = 1 + b[0]*t + b[1]*t*t + b[2]*t*t*t
        return sign * (t - num/den)

    theoretical = [norm_ppf((i + 1 - 0.375) / (n + 0.25)) for i in range(n)]

    # Standardize residuals
    std_residuals = [(r - mean_r) / std_r if std_r > 0 else 0 for r in sorted_res]

    all_x = theoretical
    all_y = std_residuals
    x_min, x_max = min(all_x) - 0.2, max(all_x) + 0.2
    y_min, y_max = min(all_y) - 0.2, max(all_y) + 0.2

    W, H = width, height
    PAD_L, PAD_R, PAD_TOP, PAD_BOT = 36, 14, 12, 28
    plot_w = W - PAD_L - PAD_R
    plot_h = H - PAD_TOP - PAD_BOT

    def px(xv):
        return PAD_L + (xv - x_min) / (x_max - x_min) * plot_w

    def py(yv):
        return PAD_TOP + (y_max - yv) / (y_max - y_min) * plot_h

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
        f'style="font-family:\'JetBrains Mono\',monospace;font-size:9px;display:block;">',
    ]

    lines.append(f'<rect x="{PAD_L}" y="{PAD_TOP}" width="{plot_w}" height="{plot_h}" fill="#F8FAFC" rx="3" stroke="#E2E8F0" stroke-width="0.5"/>')

    # Reference line y=x
    ref_min = max(x_min, y_min)
    ref_max = min(x_max, y_max)
    lines.append(
        f'<line x1="{px(ref_min):.1f}" y1="{py(ref_min):.1f}" '
        f'x2="{px(ref_max):.1f}" y2="{py(ref_max):.1f}" '
        f'stroke="#CBD5E1" stroke-width="1" stroke-dasharray="5,3"/>'
    )

    # Points
    for xv, yv in zip(all_x, all_y):
        lines.append(f'<circle cx="{px(xv):.1f}" cy="{py(yv):.1f}" r="3" fill="#1A4B8C" opacity="0.75"/>')

    # Axis ticks
    for tick in [-2, -1, 0, 1, 2]:
        if x_min <= tick <= x_max:
            tx = px(tick)
            lines.append(f'<line x1="{tx:.1f}" y1="{PAD_TOP + plot_h}" x2="{tx:.1f}" y2="{PAD_TOP + plot_h + 3}" stroke="#CBD5E1" stroke-width="1"/>')
            lines.append(f'<text x="{tx:.1f}" y="{PAD_TOP + plot_h + 13}" text-anchor="middle" fill="#94A3B8">{tick}</text>')
        if y_min <= tick <= y_max:
            ty = py(tick)
            lines.append(f'<line x1="{PAD_L - 3}" y1="{ty:.1f}" x2="{PAD_L}" y2="{ty:.1f}" stroke="#CBD5E1" stroke-width="1"/>')
            lines.append(f'<text x="{PAD_L - 5}" y="{ty + 3:.1f}" text-anchor="end" fill="#94A3B8">{tick}</text>')

    # Labels
    lines.append(f'<text x="{PAD_L + plot_w // 2}" y="{H - 2}" text-anchor="middle" fill="#64748B" font-size="9px">Theoretical Quantiles</text>')
    lines.append(f'<text x="9" y="{PAD_TOP + plot_h // 2}" text-anchor="middle" fill="#64748B" font-size="9px" transform="rotate(-90,9,{PAD_TOP + plot_h // 2})">Std Residuals</text>')

    lines.append('</svg>')
    return '\n'.join(lines)
