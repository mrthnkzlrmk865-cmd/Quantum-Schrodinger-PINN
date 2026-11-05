"""
Scattering observables for the open disordered waveguide
(config.waveguide_config): a wave packet in the lowest transverse mode
hits a disordered slab a < x < b; the CAPs at both ends absorb outgoing
waves.

At time t the probability splits exactly into
    T(t)  = absorbed_right(t) + P(x > b, t)   transmitted
    R(t)  = absorbed_left(t)  + P(x < a, t)   reflected (or not yet arrived)
    P_in(t) = P(a < x < b, t)                 still inside the slab
with T + R + P_in = 1 (up to rounding). Once the packet has left the
source region, T and R are cumulative transmission and reflection
probabilities. If P_in(t_end) is not small, the
asymptotic T and R are not yet reached; this is reported explicitly,
not hidden. A slowly decaying P_in(t) is also a physical observable
(the dwell or trapping time in the disordered region).

The observables are averaged over the packet's energy distribution
(sigma_k ~ 1/(sqrt 2 sigma_x)). They are not the fixed-energy transmission
eigenvalues of mesoscopic theory.
"""

from __future__ import annotations

import time

import numpy as np

from disorder import disorder_from_config
from potential import Potential
from reference_solver import SplitStepDST


def run_transport(cfg, seeds, t_end=None, n_frames=None, workers=-1, batch=16):
    """Propagate one realization per seed. Returns dict of arrays:
    times (F,), T/R/P_in (F, n_seeds), plus seeds and wall time."""
    t_end = cfg.T if t_end is None else t_end
    a, b = cfg.disorder_region
    seeds = list(seeds)
    out = {k: [] for k in ("T", "R", "P_in")}
    t0 = time.time()
    times = None
    for s in range(0, len(seeds), batch):
        chunk = seeds[s:s + batch]
        pots = [Potential(cfg, disorder_from_config(cfg, seed=sd)) for sd in chunk]
        solver = SplitStepDST(cfg, pots, workers=workers)
        xs = solver.X
        right = xs > b
        left = xs < a
        dxdy = solver.dx * solver.dy

        def obs(psi, slv):
            rho = np.abs(psi) ** 2
            return dict(P_right=(rho * right).sum(axis=(-2, -1)) * dxdy,
                        P_left=(rho * left).sum(axis=(-2, -1)) * dxdy)

        res = solver.run(t_end=t_end, n_frames=n_frames, store_frames=False, observables=obs)
        times = res["times"]
        T = res["absorbed_right"] + res["P_right"]
        R = res["absorbed_left"] + res["P_left"]
        out["T"].append(T)
        out["R"].append(R)
        out["P_in"].append(1.0 - T - R)
    res = {k: np.concatenate(v, axis=1) for k, v in out.items()}
    res["times"] = times
    res["seeds"] = np.array(seeds)
    res["wall_time_s"] = time.time() - t0
    return res
