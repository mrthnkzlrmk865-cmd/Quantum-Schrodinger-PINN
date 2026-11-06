"""
B1 clean-waveguide validation and numerical-stability study of the
stationary RGF solver (stationary.py). Results: docs/RGF_B1_VALIDATION.md.

B1 thresholds (docs/PHASE2_METHOD_SELECTION.md §3, unchanged):
  |T_nn - 1| <= 1e-12; sum_{m!=n} T_mn <= 1e-12; R_n <= 1e-12;
  U = max(|S^dag S - I|, |S S^dag - I|) <= 1e-10; phase |t_nn - lam^(N-1)| <= 1e-10
  (provisional); reciprocity Q = max|S - S^T| <= 1e-8 (value of B2/B3; B1 has none).

Studies:
  b1_grid       : spec B1 energy grid (50 energies, 1st-6th channel windows,
                  |E - E_perp| >= 1e-3 E) for N_y in {63, 127}, dx in {1/512, 1/1024}
                  on the B3 device (slab 1.9-2.9 + >= 0.2 buffers), and for the
                  longest L_d = 3 device of the length sweep at both dx
  near_threshold: E = E_perp_n (1 +/- delta), delta in {1e-3, 1e-5, 1e-7, 1e-9},
                  n = 1..7 (outside the B1 grid; measured, not pass/fail)
  phase_exact   : the phase error against the exact lattice factor lam^(N-1)
                  evaluated in 50-digit arithmetic (mpmath)
Writes results/rgf_b1_validation/{b1_grid.csv, near_threshold.csv,
phase_exact.csv, summary.json, run_metadata.json}.

Usage: python validate_rgf_b1.py [--out DIR] [--quick]
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")    # performance (stationary.py docstring); before NumPy

import argparse
import csv
import json
import time

import numpy as np

import smatrix_spec as sm
from config import waveguide_config
from potential import Potential
from provenance import make_metadata, write_dir_metadata
from stationary import build_device, scattering_matrix

LIMITS = dict(T_diag=1e-12, mixing=1e-12, R=1e-12, U=1e-10, phase=1e-10, Q=1e-8)


def b1_errors(r, N):
    lam = r["modes"]["lam"][r["modes"]["open"]]
    S = r["S"]
    T, Tb = np.abs(r["t_RL"]) ** 2, np.abs(r["t_LR"]) ** 2
    I = np.eye(S.shape[0])
    return dict(
        T_diag=float(max(np.max(np.abs(np.diag(T) - 1)), np.max(np.abs(np.diag(Tb) - 1)))),
        mixing=float(max(np.max(T.sum(0) - np.diag(T)), np.max(Tb.sum(0) - np.diag(Tb)))),
        R=float(max(np.max((np.abs(r["r_L"]) ** 2).sum(0)), np.max((np.abs(r["r_R"]) ** 2).sum(0)))),
        U=float(max(np.max(np.abs(S.conj().T @ S - I)), np.max(np.abs(S @ S.conj().T - I)))),
        phase=float(max(np.max(np.abs(np.diag(r["t_RL"]) - lam ** (N - 1))),
                        np.max(np.abs(np.diag(r["t_LR"]) - lam ** (N - 1))))),
        Q=float(sm.reciprocity_residual(S)))


def spec_energies(E_perp):
    Es = np.linspace(E_perp[0], E_perp[6], 52)[1:-1]
    return np.array([E for E in Es if np.all(np.abs(E - E_perp) >= 1e-3 * E)])


def device(ny, dx, slab_length):
    """Clean guide on the region of the transport device for slab length
    L_d: slab 1.9 .. 1.9 + L_d with >= 0.2 buffers, rounded up to whole dx."""
    pot = Potential(waveguide_config(slab_length=slab_length))     # V = 0 (no disorder, no barrier)
    x0 = 1.7
    n = int(np.ceil((slab_length + 0.4) / dx - 1e-9))
    return build_device(pot.numpy, x0, x0 + n * dx, dx, 1.0, ny)


def run_b1_grid(configs):
    rows = []
    for ny, dx, Ld in configs:
        dev = device(ny, dx, Ld)
        Es = spec_energies(dev.basis.E_perp)
        t0 = time.time()
        for E in Es:
            r = scattering_matrix(E, dev)
            err = b1_errors(r, dev.n_slices)
            rows.append(dict(N_y=ny, dx=dx, L_d=Ld, N=dev.n_slices, E=float(E), n_open=r["n_open"],
                             n_closed=ny - r["n_open"], min_rcond=r["min_rcond"], **err))
        print(f"  b1 grid N_y={ny} dx=1/{round(1 / dx)} L_d={Ld} N={dev.n_slices}: {len(Es)} energies "
              f"in {time.time() - t0:.1f}s", flush=True)
    return rows


def run_near_threshold(ny=63, dx=1 / 512, Ld=1.0):
    dev = device(ny, dx, Ld)
    rows = []
    for n in range(1, 8):
        En = float(dev.basis.E_perp[n - 1])
        for sgn in (-1, +1):
            for delta in (1e-3, 1e-5, 1e-7, 1e-9):
                E = En * (1 + sgn * delta)
                if n == 1 and sgn < 0:
                    continue                         # below the first threshold: no open channel
                try:
                    r = scattering_matrix(E, dev)
                except ValueError as exc:            # lead_modes rejects |c| = 1 within 1e-12
                    rows.append(dict(n=n, side=sgn, delta=delta, E=E, status=f"rejected: {exc}"))
                    continue
                m = r["modes"]
                sin_min = float(np.min(np.sin(m["k"][m["open"]] * dx)))
                rows.append(dict(n=n, side=sgn, delta=delta, E=E, status="ok", n_open=r["n_open"],
                                 min_sin_kdx=sin_min, min_rcond=r["min_rcond"], **b1_errors(r, dev.n_slices)))
    return rows


def run_phase_exact(ny=63, dx=1 / 512, Ld=1.0):
    import mpmath as mp
    mp.mp.dps = 50
    dev = device(ny, dx, Ld)
    N, Ep = dev.n_slices, dev.basis.E_perp
    rows = []
    for E in spec_energies(Ep):
        r = scattering_matrix(E, dev)
        o = r["modes"]["open"]
        for n in np.where(o)[0]:
            lam = r["modes"]["lam"][n]
            c = 1 - (mp.mpf(float(E)) - mp.mpf(float(Ep[n]))) * mp.mpf(dx) ** 2
            exact = complex(mp.exp(1j * mp.acos(c)) ** (N - 1))
            k = list(np.where(o)[0]).index(n)
            rows.append(dict(E=float(E), n=int(n) + 1, sin_kdx=float(np.sin(np.angle(lam))),
                             err_rgf_vs_exact=abs(r["t_RL"][k, k] - exact),
                             err_double_ref_vs_exact=abs(lam ** (N - 1) - exact),
                             err_rgf_vs_double_ref=abs(r["t_RL"][k, k] - lam ** (N - 1))))
    return rows


def write(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rgf_b1_validation")
    ap.add_argument("--quick", action="store_true", help="only N_y = 63, dx = 1/512 (smoke run)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    configs = [(63, 1 / 512, 1.0)] if a.quick else [
        (63, 1 / 512, 1.0), (127, 1 / 512, 1.0), (63, 1 / 1024, 1.0), (127, 1 / 1024, 1.0),
        (63, 1 / 512, 3.0), (63, 1 / 1024, 3.0)]
    grid = run_b1_grid(configs)
    near = run_near_threshold()
    phase = run_phase_exact()
    write(os.path.join(a.out, "b1_grid.csv"), grid)
    write(os.path.join(a.out, "near_threshold.csv"), near)
    write(os.path.join(a.out, "phase_exact.csv"), phase)
    per_cfg = {}
    for ny, dx, Ld in configs:
        rr = [r for r in grid if r["N_y"] == ny and r["dx"] == dx and r["L_d"] == Ld]
        mx = {k: max(r[k] for r in rr) for k in LIMITS}
        per_cfg[f"N_y={ny},dx=1/{round(1 / dx)},L_d={Ld},N={rr[0]['N']}"] = dict(
            n_energies=len(rr), max=mx, min_rcond=min(r["min_rcond"] for r in rr),
            n_open_range=[min(r["n_open"] for r in rr), max(r["n_open"] for r in rr)],
            passed={k: mx[k] <= LIMITS[k] for k in LIMITS})
    ok_near = [r for r in near if r["status"] == "ok"]
    summary = dict(
        limits=LIMITS, b1_grid=per_cfg,
        all_b1_passed=all(all(v["passed"].values()) for v in per_cfg.values()),
        near_threshold=dict(
            n_ok=len(ok_near), n_rejected=len(near) - len(ok_near),
            max_by_delta={str(d): {k: max(r[k] for r in ok_near if r["delta"] == d) for k in LIMITS}
                          for d in (1e-3, 1e-5, 1e-7, 1e-9) if any(r["delta"] == d for r in ok_near)},
            min_rcond=min(r["min_rcond"] for r in ok_near)),
        phase_exact=dict(max_err_rgf_vs_exact=max(r["err_rgf_vs_exact"] for r in phase),
                         max_err_double_ref_vs_exact=max(r["err_double_ref_vs_exact"] for r in phase)),
        openblas_num_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
        wall_time_s=time.time() - t0)
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    write_dir_metadata(a.out, make_metadata("rgf", dict(experiment="rgf_b1_validation",
                                                       configs=[list(c) for c in configs]),
                                            wall_time_s=summary["wall_time_s"],
                                            openblas_num_threads=summary["openblas_num_threads"]))
    print(json.dumps({k: summary[k] for k in ("all_b1_passed", "near_threshold", "phase_exact", "wall_time_s")},
                     indent=1))
