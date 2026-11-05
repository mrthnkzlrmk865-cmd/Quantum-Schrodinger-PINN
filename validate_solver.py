"""
Numerical validation of the DST split-step reference solver.

The reference solver is NOT assumed to be ground truth. This script measures
its accuracy for both configurations used in the project:

  V1  exact eigenmode evolution (constant V): phase error vs time
  V2  free Gaussian packet vs the analytic solution (2x2 box, walls far away)
  V3  time-step convergence (Richardson / self-convergence) with strong
      disorder: PINN problem (wavefunction) and waveguide (T, R)
  V4  spatial convergence with the SAME continuum disorder realization on
      dx = 1/32, 1/64, 1/128 (waveguide T, R; PINN-problem wavefunction)
  V5  spurious CAP return vs incident momentum (clean waveguide)
  V6  norm conservation (closed box) and exact flux bookkeeping (open)
  V7  dependence of ensemble-mean transmission on the number of disorder
      Fourier modes M

Writes results/validation/*.csv and validation_summary.json.
Usage:  python validate_solver.py            (~15-25 min on 4 CPU cores)
"""

from __future__ import annotations

import csv
import json
import os
import time

import numpy as np

from analytic import box_eigenmode, free_gaussian_2d
from config import pinn_config, waveguide_config
from potential import Potential
from reference_solver import SplitStepDST
from transport import run_transport

OUT = "results/validation"


def write_csv(name, rows):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


def v1_eigenmode():
    cfg = pinn_config(nx=63, ny=63, T=0.1)
    s = SplitStepDST(cfg, np.full((63, 63), 50.0))
    psi0 = box_eigenmode(s.X, s.Y, 0, 1, 1, 3, 2)
    res = s.run(psi0=psi0, dt=1e-4, n_frames=11)
    rows = [dict(t=t, max_abs_err=float(np.abs(f[0] - box_eigenmode(s.X, s.Y, t, 1, 1, 3, 2, V0=50.0)).max()))
            for t, f in zip(res["times"], res["frames"])]
    write_csv("v1_eigenmode.csv", rows)
    return max(r["max_abs_err"] for r in rows)


def v2_free_gaussian():
    rows = []
    for n in (127, 255):
        cfg = pinn_config(T=0.01, Lx=2.0, Ly=2.0, x0=0.85, y0=1.0, nx=n, ny=n)
        s = SplitStepDST(cfg, Potential(cfg))
        for dt in (1e-4, 2e-5):
            res = s.run(dt=dt, n_frames=6)
            for t, f in zip(res["times"], res["frames"][:, 0]):
                ex = free_gaussian_2d(cfg, s.X, s.Y, t)
                rows.append(dict(n=n, dt=dt, t=t, rel_l2=float(np.linalg.norm(f - ex) / np.linalg.norm(ex))))
    write_csv("v2_free_gaussian.csv", rows)
    return max(r["rel_l2"] for r in rows)


def v3_time_convergence():
    rows = []
    # (a) PINN problem with W = E0 (strongest disorder used for the PINN)
    cfg = pinn_config(disorder_strength=112.5, disorder_corr_length=0.08)
    s = SplitStepDST(cfg)
    ref = s.run(dt=1.25e-6, n_frames=2)["frames"][-1, 0]
    for dt in (8e-5, 4e-5, 2e-5, 1e-5, 5e-6):
        f = s.run(dt=dt, n_frames=2)["frames"][-1, 0]
        rows.append(dict(case="pinn_W112", dt=dt, quantity="rel_l2_psi(T)",
                         error=float(np.linalg.norm(f - ref) / np.linalg.norm(ref))))
    # (b) waveguide, strong disorder W = 200, one realization, observables at t = 1
    cfg = waveguide_config(disorder_strength=200.0)
    vals = {}
    for dt in (8e-4, 4e-4, 2e-4, 1e-4, 2.5e-5):
        r = run_transport(cfg.replace(dt=dt), [0], t_end=1.0, n_frames=2)
        vals[dt] = (r["T"][-1, 0], r["R"][-1, 0])
    T_ref, R_ref = vals[2.5e-5]
    for dt in (8e-4, 4e-4, 2e-4, 1e-4):
        rows.append(dict(case="waveguide_W200", dt=dt, quantity="|T-T_ref|", error=float(abs(vals[dt][0] - T_ref))))
        rows.append(dict(case="waveguide_W200", dt=dt, quantity="|R-R_ref|", error=float(abs(vals[dt][1] - R_ref))))
    write_csv("v3_time_convergence.csv", rows)
    return rows


def v4_space_convergence():
    rows = []
    # waveguide: same continuum potential sampled on three grids
    base = waveguide_config(disorder_strength=200.0)
    for zeta in (0.02, 0.05):
        for npl in (32, 64, 128):
            nx = int(round(base.Lx * npl)) - 1
            cfg = base.replace(nx=nx, ny=npl - 1, disorder_corr_length=zeta, dt=1e-4)
            r = run_transport(cfg, [0, 1], t_end=1.0, n_frames=2)
            for i, sd in enumerate((0, 1)):
                rows.append(dict(case=f"waveguide_W200_zeta{zeta}", dx=1.0 / npl, seed=sd,
                                 T=float(r["T"][-1, i]), R=float(r["R"][-1, i]), P_in=float(r["P_in"][-1, i])))
    # PINN problem: psi(T) at the common coarse grid points
    cfg0 = pinn_config(disorder_strength=112.5)
    fields = {}
    for n in (31, 63, 127, 255):
        cfg = cfg0.replace(nx=n, ny=n)
        f = SplitStepDST(cfg).run(dt=5e-6, n_frames=2)["frames"][-1, 0]
        fields[n] = f
    fine = fields[255]
    for n in (31, 63, 127):
        stride = 256 // (n + 1)
        sub = fine[stride - 1::stride, stride - 1::stride]
        rows.append(dict(case="pinn_W112_psi", dx=1.0 / (n + 1), seed=0,
                         T=float(np.linalg.norm(fields[n] - sub) / np.linalg.norm(sub)), R=np.nan, P_in=np.nan))
    write_csv("v4_space_convergence.csv", rows)
    return rows


def v5_cap_reflection():
    rows = []
    for k in (10.0, 12.0, 16.0, 20.0, 25.0, 30.0):
        cfg = waveguide_config(kx0=k)
        # long enough for the slowest relevant components (k - 3 sigma_k) to be absorbed
        t_end = 2.6 / max(k - 3 * 3.54, 3.0) + 0.3
        r = run_transport(cfg, [0], t_end=t_end, n_frames=2)
        rows.append(dict(kx0=k, t_end=t_end, spurious_return=float(r["R"][-1, 0]),
                         T=float(r["T"][-1, 0]), remaining=float(r["P_in"][-1, 0])))
    write_csv("v5_cap_reflection.csv", rows)
    return rows


def v6_conservation():
    cfg = pinn_config(disorder_strength=112.5)
    r = SplitStepDST(cfg).run(n_frames=21, store_frames=False)
    closed = float(np.max(np.abs(r["norm"] - 1.0)))
    cfg = waveguide_config(disorder_strength=200.0)
    s = SplitStepDST(cfg)
    r = s.run(t_end=1.0, n_frames=11, store_frames=False)
    book = float(np.max(np.abs(r["norm"] + r["absorbed_left"] + r["absorbed_right"] - r["norm"][0])))
    return dict(closed_norm_max_dev=closed, open_bookkeeping_max_dev=book)


def v7_n_modes(n_real=24):
    rows = []
    for M in (64, 256, 512, 2048):
        cfg = waveguide_config(disorder_strength=100.0, disorder_n_modes=M)
        r = run_transport(cfg, range(n_real), t_end=1.0, n_frames=2)
        T = r["T"][-1]
        rows.append(dict(n_modes=M, n_real=n_real, T_mean=float(T.mean()),
                         T_se=float(T.std(ddof=1) / np.sqrt(n_real))))
    write_csv("v7_n_modes.csv", rows)
    return rows


if __name__ == "__main__":
    t0 = time.time()
    summary = {}
    summary["v1_eigenmode_max_abs_err"] = v1_eigenmode(); print("V1", summary["v1_eigenmode_max_abs_err"], flush=True)
    summary["v2_free_gaussian_max_rel_l2"] = v2_free_gaussian(); print("V2", summary["v2_free_gaussian_max_rel_l2"], flush=True)
    summary["v3"] = v3_time_convergence(); print("V3", summary["v3"], flush=True)
    summary["v4"] = v4_space_convergence(); print("V4", summary["v4"], flush=True)
    summary["v5"] = v5_cap_reflection(); print("V5", summary["v5"], flush=True)
    summary["v6"] = v6_conservation(); print("V6", summary["v6"], flush=True)
    summary["v7"] = v7_n_modes(); print("V7", summary["v7"], flush=True)
    summary["wall_time_s"] = time.time() - t0
    with open(os.path.join(OUT, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"done in {summary['wall_time_s']:.0f}s")
