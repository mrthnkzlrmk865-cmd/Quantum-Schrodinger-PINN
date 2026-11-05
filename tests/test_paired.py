"""Paired-bootstrap analysis of the correlation-length sweep."""
import json
import os

import numpy as np
import pytest

import paired_analysis as pa


def test_paired_contrast_mean_and_pairing_benefit():
    rng = np.random.default_rng(0)
    common = rng.normal(0.5, 0.1, 48)                 # strongly shared per-seed component
    a = common + rng.normal(0, 0.01, 48)
    b = common + 0.02 + rng.normal(0, 0.01, 48)
    idx = pa.paired_boot_indices(48, n_boot=4000)
    c = pa.paired_contrast({"a": a, "b": b}, {"a": 1.0, "b": -1.0}, idx)
    assert np.isclose(c["mean"], (a - b).mean())
    paired = np.diff(c["boot_ci_0.9500"])[0]
    unpaired = np.diff(c["unpaired_welch_ci_0.95"])[0]
    assert paired < unpaired / 5                      # pairing removes the shared variance
    assert c["r_pairs"] > 0.9
    assert c["boot_ci_0.9833"][0] <= c["boot_ci_0.9500"][0]   # Bonferroni interval is wider


def test_paired_bootstrap_coverage():
    rng = np.random.default_rng(1)
    hits = 0
    for _ in range(300):
        d = rng.normal(0.03, 0.05, 48)
        idx = pa.paired_boot_indices(48, n_boot=1000, seed=int(rng.integers(1 << 30)))
        lo, hi = pa.paired_contrast({"x": d, "z": np.zeros(48)}, {"x": 1.0, "z": -1.0}, idx)["boot_ci_0.9500"]
        hits += lo <= 0.03 <= hi
    assert 0.90 <= hits / 300 <= 0.98


def test_curvature_contrast_and_vertex():
    z = np.array([0.03, 0.05, 0.07])
    lx = np.log(z)
    v_true = np.log(0.045)
    T = {zz: np.full(10, 0.5 + 2.0 * (np.log(zz) - v_true) ** 2) for zz in z}
    idx = pa.paired_boot_indices(10, n_boot=200)
    K = pa.paired_contrast(T, {0.03: 0.5, 0.07: 0.5, 0.05: -1.0}, idx)["mean"]
    assert K > 0
    assert np.isclose(np.exp(pa.parabola_vertex(lx, [T[zz][0] for zz in z])), 0.045)
    assert np.isnan(pa.parabola_vertex(lx, [0.5, 0.6, 0.5]))   # concave: no minimum


def test_holm():
    adj = pa.holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert np.isclose(adj["a"], 0.03) and np.isclose(adj["c"], 0.06) and np.isclose(adj["b"], 0.06)


def _fake_point(tmp, z, seeds, cfg_extra=None):
    cfg = dict(disorder_corr_length=z, kx0=20.0, Ly=1.0, disorder_strength=100.0)
    cfg.update(cfg_extra or {})
    meta = dict(spec=dict(config=cfg, run=dict(t_end=1.0, n_real=len(seeds))))
    np.savez(os.path.join(tmp, f"zeta{z}.npz"), T=np.ones((2, len(seeds))), seeds=np.asarray(seeds),
             metadata_json=json.dumps(meta))


def test_load_sweep_refuses_unpaired_or_inconsistent_data(tmp_path):
    _fake_point(tmp_path, 0.02, [0, 1, 2])
    _fake_point(tmp_path, 0.05, [0, 1, 3])
    with pytest.raises(ValueError, match="seeds differ"):
        pa.load_sweep(str(tmp_path))
    for f in tmp_path.iterdir():
        f.unlink()
    _fake_point(tmp_path, 0.02, [0, 1, 2])
    _fake_point(tmp_path, 0.05, [0, 1, 2], dict(disorder_strength=101.0))
    with pytest.raises(ValueError, match="disorder_strength"):
        pa.load_sweep(str(tmp_path))


def test_load_sweep_on_reproduced_data():
    T, seeds, k = pa.load_sweep("results/reproduction/transport/corr_length")
    assert sorted(T) == [0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2]
    assert len(seeds) == 48 and np.isclose(k, np.hypot(20, np.pi))
    assert np.isclose(T[0.05].mean(), 0.494, atol=1e-3)


# ----------------------------------------------------------------------
# bootstrap p-value with the finite-resampling (+1) correction
# ----------------------------------------------------------------------
def test_bootstrap_p_value_plus_one_correction():
    boot = np.full(20000, 0.10) + np.linspace(-0.01, 0.01, 20000)   # all near the observation
    r = pa.bootstrap_p_value(boot, observed=0.10)
    assert r["p_boot_count"] == 0
    assert r["p_boot"] == 1 / 20001 and r["p_boot"] > 0                # never exactly zero
    assert r["p_boot_resolution"] == 1 / 20001 and r["p_boot_B"] == 20000


def test_bootstrap_p_value_counting_and_centering():
    obs = 0.5
    boot = np.array([0.5, 0.0, 1.0, 1.2, 0.1, 0.9, -0.2])
    # centered replicates boot - obs; extreme if |boot - obs| >= |obs| = 0.5
    expected_b = int(np.sum(np.abs(boot - obs) >= 0.5))
    assert expected_b == 4                                       # 0.0, 1.0, 1.2, -0.2
    r = pa.bootstrap_p_value(boot, obs)
    assert r["p_boot_count"] == 4 and np.isclose(r["p_boot"], 5 / 8)


def test_bootstrap_p_value_null_observation_is_one():
    r = pa.bootstrap_p_value(np.random.default_rng(0).normal(size=1000), 0.0)
    assert r["p_boot"] == 1.0                                   # every replicate is "as extreme"


def test_bootstrap_p_value_finite_resolution_and_holm_floor():
    """With B replicates no p below 1/(B+1) can be reported; after Holm over
    three contrasts the smallest reportable value is 3/(B+1)."""
    B = 20000
    pmin = pa.bootstrap_p_value(np.zeros(B) + 1.0, 1.0)["p_boot"]
    adj = pa.holm({"C1": pmin, "C2": pmin, "C3": 0.01})
    assert np.isclose(adj["C1"], 3 / (B + 1)) and np.isclose(adj["C2"], 3 / (B + 1))


def test_paired_contrast_reports_plus_one_fields():
    rng = np.random.default_rng(2)
    a = rng.normal(0.5, 0.05, 48)
    idx = pa.paired_boot_indices(48, n_boot=500)
    c = pa.paired_contrast({"a": a, "b": a - 0.2}, {"a": 1.0, "b": -1.0}, idx)
    assert c["p_boot"] == (c["p_boot_count"] + 1) / (500 + 1)
    assert c["p_boot"] >= 1 / 501


def test_analysis_settings_preserved():
    """Seed set, contrasts, bootstrap seed, B and the primary family are fixed."""
    assert pa.N_BOOT == 20000 and pa.BOOT_SEED == 20261105
    assert pa.PRIMARY == ("C1", "C2", "C3")


def test_analysis_deterministic_and_plus_one_on_real_data():
    res1, _ = pa.analyse("results/reproduction/transport/corr_length")
    res2, _ = pa.analyse("results/reproduction/transport/corr_length")
    assert res1 == res2
    for name in pa.PRIMARY:
        c = res1["contrasts"][name]
        assert c["p_boot"] == (c["p_boot_count"] + 1) / (pa.N_BOOT + 1) > 0
        assert "p_boot_holm" in c
