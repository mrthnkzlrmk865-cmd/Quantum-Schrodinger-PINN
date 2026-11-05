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
