import numpy as np

from stats_utils import running_mean, summarize


def test_summary_coverage():
    rng = np.random.default_rng(1)
    hits = 0
    for _ in range(400):
        s = summarize(rng.normal(2.0, 1.0, 30), n_boot=500)
        hits += s["ci_t_lo"] < 2.0 < s["ci_t_hi"]
    assert 0.91 < hits / 400 < 0.98


def test_running_mean_final_value():
    v = np.arange(10.0)
    n, m, se = running_mean(v)
    assert m[-1] == v.mean() and np.isclose(se[-1], v.std(ddof=1) / np.sqrt(10))
