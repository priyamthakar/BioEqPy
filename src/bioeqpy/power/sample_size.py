"""Sample-size and power estimation for 2x2 crossover BE studies."""

from __future__ import annotations

import math
from collections.abc import Iterable

import pandas as pd
from scipy import stats


def _sigma_from_cv(cv: float) -> float:
    if cv <= 0:
        raise ValueError("CV must be positive.")
    return math.sqrt(math.log(cv**2 + 1.0))


def achieved_power(
    n: int,
    cv: float,
    gmr: float = 0.95,
    alpha: float = 0.05,
    theta: tuple[float, float] = (0.80, 1.25),
) -> float:
    """Approximate TOST power for a standard balanced 2x2 crossover study."""
    if n < 4:
        raise ValueError("2x2 crossover sample size must be at least 4.")
    sigma = _sigma_from_cv(cv)
    df = n - 2
    se = sigma * math.sqrt(2.0 / n)
    delta = math.log(gmr)
    tcrit = stats.t.ppf(1.0 - alpha, df)
    lower_ncp = (delta - math.log(theta[0])) / se
    upper_ncp = (math.log(theta[1]) - delta) / se
    lower_power = 1.0 - stats.nct.cdf(tcrit, df, lower_ncp)
    upper_power = 1.0 - stats.nct.cdf(tcrit, df, upper_ncp)
    return float(max(0.0, min(1.0, lower_power + upper_power - 1.0)))


def sample_size(
    design: str = "2x2",
    cv: float = 0.25,
    gmr: float = 0.95,
    power: float = 0.80,
    alpha: float = 0.05,
    dropout_rate: float = 0.0,
    max_n: int = 1000,
) -> int:
    """Return the smallest even N meeting the target approximate TOST power."""
    if design != "2x2":
        raise NotImplementedError("Only 2x2 crossover sample size is implemented in v0.1.")
    if not 0 <= dropout_rate < 1:
        raise ValueError("Dropout rate must be in [0, 1).")
    for n in range(4, max_n + 1, 2):
        if achieved_power(n=n, cv=cv, gmr=gmr, alpha=alpha) >= power:
            adjusted = math.ceil(n / (1.0 - dropout_rate))
            return adjusted if adjusted % 2 == 0 else adjusted + 1
    raise ValueError(f"Target power was not reached by N={max_n}.")


def power_curve(
    n_range: Iterable[int],
    cv: float = 0.25,
    gmr: float = 0.95,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Build a power curve table for candidate sample sizes."""
    return pd.DataFrame(
        [{"N": int(n), "Power": achieved_power(int(n), cv=cv, gmr=gmr, alpha=alpha)} for n in n_range]
    )

