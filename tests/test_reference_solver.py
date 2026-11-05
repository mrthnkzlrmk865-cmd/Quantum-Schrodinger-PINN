"""Reference solver: exactness, convergence, conservation, CAP bookkeeping."""
import numpy as np

from analytic import box_eigenmode, free_gaussian_2d
from config import pinn_config, waveguide_config
from reference_solver import SplitStepDST
from potential import Potential


def test_eigenmode_phase_is_exact():
    cfg = pinn_config(Lx=1.0, Ly=1.3, nx=31, ny=41, T=0.05)
    s = SplitStepDST(cfg, np.full((cfg.nx, cfg.ny), 7.0))  # constant V
    psi0 = box_eigenmode(s.X, s.Y, 0.0, cfg.Lx, cfg.Ly, 2, 3)
    res = s.run(psi0=psi0, dt=1e-3, n_frames=2)
    exact = box_eigenmode(s.X, s.Y, cfg.T, cfg.Lx, cfg.Ly, 2, 3, V0=7.0)
    assert np.max(np.abs(res["frames"][-1, 0] - exact)) < 1e-10


def test_norm_conserved_without_cap():
    cfg = pinn_config(disorder_strength=112.0, nx=63, ny=63)
    res = SplitStepDST(cfg).run(dt=1e-4, n_frames=5)
    assert np.max(np.abs(res["norm"] - res["norm"][0])) < 1e-12


def test_free_gaussian_matches_analytic():
    # 2x2 box with the packet in the centre: the walls stay > 6 widths away
    # up to t = 0.01, so the unbounded analytic solution applies.
    cfg = pinn_config(T=0.01, Lx=2.0, Ly=2.0, x0=0.85, y0=1.0, nx=255, ny=255)
    s = SplitStepDST(cfg, Potential(cfg))
    res = s.run(dt=2e-5, n_frames=3)
    for t, f in zip(res["times"], res["frames"][:, 0]):
        ex = free_gaussian_2d(cfg, s.X, s.Y, t)
        assert np.linalg.norm(f - ex) / np.linalg.norm(ex) < 1e-8


def test_second_order_in_time_with_disorder():
    cfg = pinn_config(disorder_strength=112.0, nx=63, ny=63, T=0.01)
    s = SplitStepDST(cfg)
    ref = s.run(dt=1.25e-6, n_frames=2)["frames"][-1, 0]
    errs = [np.linalg.norm(s.run(dt=dt, n_frames=2)["frames"][-1, 0] - ref) for dt in (4e-5, 2e-5, 1e-5)]
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(orders > 1.8) and np.all(orders < 2.2), orders


def test_cap_bookkeeping_and_low_reflection():
    cfg = waveguide_config()
    res = SplitStepDST(cfg).run(t_end=0.4, n_frames=5, store_frames=False)
    total = res["norm"] + res["absorbed_left"] + res["absorbed_right"]
    assert np.max(np.abs(total - 1.0)) < 1e-10
    # free waveguide: all probability ends in the right CAP
    assert res["absorbed_left"][-1, 0] < 1e-3
    assert res["absorbed_right"][-1, 0] > 0.99


def test_batched_equals_individual():
    cfg = pinn_config(disorder_strength=50.0, nx=31, ny=31, T=0.005)
    from disorder import disorder_from_config
    pots = [Potential(cfg, disorder_from_config(cfg, seed=s)) for s in (0, 1)]
    b = SplitStepDST(cfg, pots).run(dt=1e-4, n_frames=2)["frames"][-1]
    for i, p in enumerate(pots):
        single = SplitStepDST(cfg, p).run(dt=1e-4, n_frames=2)["frames"][-1, 0]
        assert np.allclose(b[i], single, atol=1e-13)
