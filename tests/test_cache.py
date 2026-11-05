"""Cache correctness: results are reused only for an identical
configuration, and a mismatching or unverifiable cache is never accepted
silently or overwritten without --recompute."""
import os

import numpy as np
import pytest

import experiment_disorder as ed
import experiment_localization as el
import experiment_pinn as ep
import provenance
from config import pinn_config, waveguide_config
from provenance import CacheMismatchError, fingerprint

# A tiny waveguide (coarse grid, a few time steps) so real runs take < 1 s.
TINY = dict(nx=63, ny=15, dt=2e-3)
T_END = 0.02


# ----------------------------------------------------------------------
# Fingerprint
# ----------------------------------------------------------------------
def test_fingerprint_canonical():
    a = dict(x=1, y=dict(b=2.5, a=[1, 2]))
    b = dict(y=dict(a=[1.0, 2.0], b=2.5), x=1.0)        # other order, int vs float
    assert fingerprint(a) == fingerprint(b)
    assert fingerprint(dict(x=(0.1, 0.2))) == fingerprint(dict(x=[0.1, 0.2]))
    assert fingerprint(dict(x=0.1)) != fingerprint(dict(x=0.1 + 1e-15))


@pytest.mark.parametrize("field,value", [
    ("dt", 1e-4), ("nx", 127), ("ny", 31), ("disorder_strength", 101.0),
    ("disorder_corr_length", 0.051), ("disorder_seed", 1), ("disorder_n_modes", 256),
    ("cap_strength", 299.0), ("kx0", 20.5), ("sigma_x", 0.21), ("barrier_height", 1.0),
    ("disorder_region", (1.9, 2.95)), ("Ly", 1.01),
])
def test_transport_spec_sensitive_to_every_physical_field(field, value):
    cfg = waveguide_config(disorder_strength=100.0)
    base = fingerprint(ed.transport_spec(cfg, 8, 1.0))
    assert fingerprint(ed.transport_spec(cfg.replace(**{field: value}), 8, 1.0)) != base


def test_transport_spec_sensitive_to_run_parameters_not_to_pinn_fields():
    cfg = waveguide_config(disorder_strength=100.0)
    base = fingerprint(ed.transport_spec(cfg, 8, 1.0))
    assert fingerprint(ed.transport_spec(cfg, 9, 1.0)) != base           # realizations / seeds
    assert fingerprint(ed.transport_spec(cfg, 8, 1.5)) != base           # propagation time
    assert fingerprint(ed.transport_spec(cfg, 8, 1.0, n_frames=21)) != base
    # training settings do not affect a reference-solver result
    assert fingerprint(ed.transport_spec(cfg.replace(epochs_adam=1, seed=7), 8, 1.0)) == base


@pytest.mark.parametrize("field,value", [
    ("epochs_adam", 1999), ("seed", 43), ("hidden_units", 65), ("use_fourier", False),
    ("sampling", "rad"), ("T", 0.021), ("disorder_seed", 3), ("lr_adam", 2e-3), ("n_interior", 4001),
])
def test_pinn_spec_sensitive(field, value):
    cfg = pinn_config(disorder_strength=56.25)
    assert fingerprint(ep.pinn_spec(cfg.replace(**{field: value}))) != fingerprint(ep.pinn_spec(cfg))


def test_pinn_spec_ignores_device_and_logging():
    cfg = pinn_config()
    assert fingerprint(ep.pinn_spec(cfg.replace(device="cuda", log_every=7))) == fingerprint(ep.pinn_spec(cfg))


# ----------------------------------------------------------------------
# Transport cache (real tiny runs)
# ----------------------------------------------------------------------
@pytest.fixture
def counted_transport(monkeypatch):
    calls = []
    real = ed.run_transport

    def spy(*a, **k):
        calls.append(1)
        return real(*a, **k)

    monkeypatch.setattr(ed, "run_transport", spy)
    return calls


def _point(tmp, W=100.0, n_real=2, t_end=T_END, **kw):
    return ed.run_point("s", "p", dict(disorder_strength=W, **TINY), n_real, t_end, out=str(tmp), **kw)


def test_transport_cache_reused_for_identical_config(tmp_path, counted_transport):
    _, d1 = _point(tmp_path)
    _, d2 = _point(tmp_path)
    assert len(counted_transport) == 1                    # second call did not recompute
    for k in ("T", "R", "P_in", "times", "seeds"):
        assert np.array_equal(d1[k], d2[k])
    meta = provenance.read_npz_metadata(tmp_path / "s" / "p.npz")
    assert meta["fingerprint"] == str(d1["fingerprint"])


@pytest.mark.parametrize("change", [dict(W=150.0), dict(n_real=3), dict(n_real=1), dict(t_end=0.03)])
def test_transport_cache_rejected_after_config_change(tmp_path, counted_transport, change):
    _point(tmp_path)
    before = (tmp_path / "s" / "p.npz").read_bytes()
    with pytest.raises(CacheMismatchError):
        _point(tmp_path, **change)
    assert (tmp_path / "s" / "p.npz").read_bytes() == before   # not overwritten
    assert len(counted_transport) == 1


def test_transport_cache_zeta_change_rejected(tmp_path, counted_transport):
    ed.run_point("s", "p", dict(disorder_strength=100.0, disorder_corr_length=0.05, **TINY), 2, T_END,
                 out=str(tmp_path))
    with pytest.raises(CacheMismatchError, match="disorder_corr_length"):
        ed.run_point("s", "p", dict(disorder_strength=100.0, disorder_corr_length=0.07, **TINY), 2, T_END,
                     out=str(tmp_path))


def test_transport_recompute_overwrites(tmp_path, counted_transport):
    _point(tmp_path)
    cfg, d = _point(tmp_path, W=150.0, recompute=True)
    assert len(counted_transport) == 2
    meta = provenance.read_npz_metadata(tmp_path / "s" / "p.npz")
    assert meta["spec"]["config"]["disorder_strength"] == 150.0
    _point(tmp_path, W=150.0)                                  # now a hit again
    assert len(counted_transport) == 2


def test_legacy_cache_without_metadata_rejected(tmp_path, counted_transport):
    os.makedirs(tmp_path / "s")
    np.savez_compressed(tmp_path / "s" / "p.npz", T=np.zeros((2, 2)), R=np.zeros((2, 2)),
                        P_in=np.zeros((2, 2)), times=np.zeros(2), seeds=np.arange(2))
    with pytest.raises(CacheMismatchError, match="no provenance metadata"):
        _point(tmp_path)
    assert len(counted_transport) == 0


def test_code_change_warns_or_fails_in_strict_mode(tmp_path, counted_transport, monkeypatch):
    _point(tmp_path)
    real = provenance.code_fingerprint
    monkeypatch.setattr(provenance, "code_fingerprint",
                        lambda kind: dict(real(kind), code_fingerprint="sha256:changed"))
    with pytest.warns(UserWarning, match="source code"):
        _point(tmp_path)
    with pytest.raises(CacheMismatchError, match="strict_code"):
        _point(tmp_path, strict_code=True)
    assert len(counted_transport) == 1


def test_transport_metadata_complete(tmp_path):
    _point(tmp_path)
    meta = provenance.read_npz_metadata(tmp_path / "s" / "p.npz")
    run = meta["spec"]["run"]
    assert run["seeds"] == [0.0, 1.0] and run["n_real"] == 2.0 and run["t_end"] == T_END
    cfg = meta["spec"]["config"]
    for k in ("nx", "ny", "dt", "disorder_strength", "disorder_corr_length", "Lx", "cap_strength"):
        assert k in cfg
    assert meta["created_utc"] and meta["code_fingerprint"].startswith("sha256:")
    assert "reference_solver.py" in meta["code_files"]
    env = meta["environment"]
    assert "git_commit" in env and "numpy" in env and "python" in env


def test_sweep_writes_summary_and_dir_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(ed, "sweep_points",
                        lambda sweep, n: [("a", dict(disorder_strength=50.0, **TINY), 2, T_END),
                                          ("b", dict(disorder_strength=0.0, **TINY), 1, T_END)])
    rows = ed.run_sweep("strength", 2, out=str(tmp_path))
    assert [r["T_n"] for r in rows] == [2, 1]
    assert all(r["fingerprint"].startswith("sha256:") for r in rows)
    meta = provenance.read_dir_metadata(tmp_path / "strength")
    assert [p["point"] for p in meta["spec"]["points"]] == ["a", "b"]


# ----------------------------------------------------------------------
# Localization cache (real tiny runs)
# ----------------------------------------------------------------------
def test_localization_cache(tmp_path, monkeypatch):
    calls = []
    real = el.run_point
    monkeypatch.setattr(el, "run_point", lambda *a, **k: calls.append(1) or real(*a, **k))
    r1 = el.run_or_load(1.0, 100, 2, 0.004, out=str(tmp_path))
    r2 = el.run_or_load(1.0, 100, 2, 0.004, out=str(tmp_path))
    assert len(calls) == 1 and np.array_equal(r1[3], r2[3]) and r1[4] == r2[4]
    with pytest.raises(CacheMismatchError):
        el.run_or_load(1.0, 100, 3, 0.004, out=str(tmp_path))        # more realizations
    with pytest.raises(CacheMismatchError):
        el.run_or_load(1.0, 100, 2, 0.006, out=str(tmp_path))        # longer propagation
    assert len(calls) == 1


# ----------------------------------------------------------------------
# PINN run cache (training replaced by a stub: only the cache logic is tested)
# ----------------------------------------------------------------------
@pytest.fixture
def stub_pinn(monkeypatch):
    calls = []
    from model import SchrodingerPINN

    def fake_train(cfg, verbose=False):
        calls.append(cfg)
        return SchrodingerPINN(cfg), [dict(loss=1.0)], dict(train_time_s=0.0, final_loss=1.0)

    class FakeSolver:
        def __init__(self, *a, **k):
            pass

        def run(self, **k):
            return {}

    monkeypatch.setattr(ep, "train", fake_train)
    monkeypatch.setattr(ep, "SplitStepDST", FakeSolver)
    monkeypatch.setattr(ep, "evaluate_run", lambda d: dict(rel_l2_final=0.5, rel_l2_mean=0.4,
                                                          rel_l2_rho_final=0.6, norm_max_dev=0.01))
    return calls


def test_pinn_cache(tmp_path, stub_pinn):
    out = str(tmp_path)
    ep.run_one("r", dict(disorder_strength=56.25), out=out)
    ep.run_one("r", dict(disorder_strength=56.25), out=out)
    assert len(stub_pinn) == 1
    with pytest.raises(CacheMismatchError):
        ep.run_one("r", dict(disorder_strength=56.25, epochs_adam=10), out=out)
    with pytest.raises(CacheMismatchError):
        ep.run_one("r", dict(disorder_strength=28.0), out=out)
    # quick runs live in their own directory and never collide with the full run
    ep.run_one("r", dict(disorder_strength=56.25), quick=True, out=out)
    assert os.path.exists(tmp_path / "r_quick" / "run_metadata.json")
    assert len(stub_pinn) == 2
    meta = provenance.read_dir_metadata(tmp_path / "r")
    assert meta["spec"]["config"]["epochs_adam"] == 2000.0


def test_pinn_legacy_run_rejected(tmp_path, stub_pinn):
    from model import SchrodingerPINN
    from train import save_run
    cfg = pinn_config()
    save_run(str(tmp_path / "r"), cfg, SchrodingerPINN(cfg), [dict(loss=1.0)], dict(train_time_s=0.0))
    with pytest.raises(CacheMismatchError, match="no provenance metadata"):
        ep.run_one("r", dict(), out=str(tmp_path))
    assert len(stub_pinn) == 0
