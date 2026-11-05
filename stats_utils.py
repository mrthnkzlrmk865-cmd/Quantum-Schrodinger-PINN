"""
Statistics over independent disorder realizations.

For each configuration we report mean, sample standard deviation (ddof=1),
standard error SE = std/sqrt(n), a Student-t 95 % confidence interval for
the mean, and a percentile bootstrap 95 % CI (fixed RNG seed so the
numbers are reproducible). The bootstrap CI is included because some
observables (e.g. transmission close to 0 or 1, or ln T) are skewed.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def summarize(values, confidence=0.95, n_boot=4000, seed=12345):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n = v.size
    mean = v.mean() if n else np.nan
    std = v.std(ddof=1) if n > 1 else np.nan
    se = std / np.sqrt(n) if n > 1 else np.nan
    if n > 1:
        tq = stats.t.ppf(0.5 + confidence / 2, n - 1)
        t_lo, t_hi = mean - tq * se, mean + tq * se
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, n, size=(n_boot, n))
        bm = v[idx].mean(axis=1)
        a = (1 - confidence) / 2
        b_lo, b_hi = np.quantile(bm, [a, 1 - a])
    else:
        t_lo = t_hi = b_lo = b_hi = np.nan
    return dict(n=n, mean=mean, std=std, se=se, ci_t_lo=t_lo, ci_t_hi=t_hi,
                ci_boot_lo=b_lo, ci_boot_hi=b_hi,
                median=np.median(v) if n else np.nan,
                min=v.min() if n else np.nan, max=v.max() if n else np.nan)


def running_mean(values):
    """Running mean and running SE after n = 1..N realizations, used to
    check convergence of ensemble averages in the number of realizations."""
    v = np.asarray(values, dtype=float)
    n = np.arange(1, v.size + 1)
    m = np.cumsum(v) / n
    m2 = np.cumsum(v**2) / n
    var = np.where(n > 1, (m2 - m**2) * n / np.maximum(n - 1, 1), np.nan)
    return n, m, np.sqrt(np.maximum(var, 0) / n)


def welch_test(a, b):
    """Two-sample Welch t-test (unequal variances). Returns (t, p)."""
    r = stats.ttest_ind(np.asarray(a), np.asarray(b), equal_var=False)
    return float(r.statistic), float(r.pvalue)
