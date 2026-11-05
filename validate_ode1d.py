"""
Measured validation of the independent 1D ODE scattering reference
(ode1d.py); results are reported in docs/ODE1D_VALIDATION.md.

Cases (B2 of docs/PHASE2_METHOD_SELECTION.md):
  sharp rectangular barriers (V0, w) in {(150, 0.1), (250, 0.1), (250, 0.3)}:
      60 energies across the B4 packet spectrum (eps = kx^2/2,
      kx in [4.09, 35.91]) plus 5 near-threshold energies; errors vs the
      closed form and vs exact plane-wave matching (amplitudes); flux error;
      lead-range extension
  smooth barrier = Potential.barrier_numpy of waveguide_config
      (V_b in {150, 250}, width 0.1, edge 0.01): rtol 1e-10 -> 1e-12 change,
      clean-range doubling change, flux error
  convergence in rtol (1e-6 ... 1e-13) at a few energies
Writes results/ode1d_validation/{sharp.csv, smooth.csv, rtol_convergence.csv,
summary.json, run_metadata.json}.
Usage: python validate_ode1d.py [--out DIR]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time

import numpy as np

from analytic import rectangular_barrier_matching, rectangular_barrier_transmission
from config import waveguide_config
from ode1d import (DEFAULT_ATOL, DEFAULT_RTOL, clean_interval, potential_barrier_callable,
                   rectangular_barrier, shift_reference, solve_scattering)
from provenance import make_metadata, write_dir_metadata

EPS_PACKET = np.linspace(0.5 * 4.09**2, 0.5 * 35.91**2, 60)


def near(V0):
    return [V0 * (1 - 1e-3), V0 * (1 - 1e-6), V0, V0 * (1 + 1e-6), V0 * (1 + 1e-3)]


def write(path, rows):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


def sharp_cases():
    rows = []
    for V0, w in [(150.0, 0.1), (250.0, 0.1), (250.0, 0.3)]:
        V, bp = rectangular_barrier(V0, 0.0, w)
        for eps in np.concatenate([EPS_PACKET, near(V0)]):
            res = solve_scattering(eps, V, 0.0, w, bp)
            Tf = rectangular_barrier_transmission(np.array([eps]), V0, w)[0]
            tm, rm = rectangular_barrier_matching(eps, V0, w)
            wide = shift_reference(solve_scattering(eps, V, -0.5, w + 0.5, bp), 0.0, w)
            rows.append(dict(V0=V0, w=w, eps=eps, T_ode=res["T"], T_formula=Tf,
                             err_T_formula=abs(res["T"] - Tf), err_t_matching=abs(res["t"] - tm),
                             err_r_matching=abs(res["r"] - rm), flux_error=res["flux_error"],
                             range_change_T=abs(wide["T"] - res["T"]),
                             range_change_t=abs(wide["t"] - res["t"]), nfev=res["nfev"]))
    return rows


def smooth_cases():
    rows = []
    for Vb in (150.0, 250.0):
        cfg = waveguide_config(barrier_height=Vb)
        V = potential_barrier_callable(cfg)
        x1 = cfg.barrier_center - cfg.barrier_width / 2
        x2 = cfg.barrier_center + cfg.barrier_width / 2
        a, b = clean_interval(V, x1, x2, tol=1e-12)
        for eps in np.concatenate([EPS_PACKET, near(Vb)]):
            base = solve_scattering(eps, V, a, b)
            loose = solve_scattering(eps, V, a, b, rtol=1e-10, atol=1e-10)
            wide = solve_scattering(eps, V, a - (x1 - a), b + (b - x2))
            rows.append(dict(V_b=Vb, edge=cfg.barrier_edge, width=cfg.barrier_width, x_left=a, x_right=b,
                             eps=eps, T_ode=base["T"], rtol_change_T=abs(loose["T"] - base["T"]),
                             range_change_T=abs(wide["T"] - base["T"]), flux_error=base["flux_error"],
                             nfev=base["nfev"]))
    return rows


def rtol_convergence():
    rows = []
    V, bp = rectangular_barrier(250.0, 0.0, 0.1)
    cfg = waveguide_config(barrier_height=250.0)
    Vs = potential_barrier_callable(cfg)
    x1 = cfg.barrier_center - cfg.barrier_width / 2
    x2 = cfg.barrier_center + cfg.barrier_width / 2
    a, b = clean_interval(Vs, x1, x2, tol=1e-12)
    for eps in (60.0, 249.0, 400.0):
        Tf = rectangular_barrier_transmission(np.array([eps]), 250.0, 0.1)[0]
        ref_s = solve_scattering(eps, Vs, a, b, rtol=1e-13, atol=1e-13)["T"]
        for tol in (1e-6, 1e-8, 1e-10, 1e-12, 1e-13):
            rs = solve_scattering(eps, V, 0.0, 0.1, bp, rtol=tol, atol=tol)
            rsm = solve_scattering(eps, Vs, a, b, rtol=tol, atol=tol)
            rows.append(dict(eps=eps, rtol=tol, sharp_err_vs_formula=abs(rs["T"] - Tf),
                             sharp_flux_error=rs["flux_error"], sharp_nfev=rs["nfev"],
                             smooth_change_vs_tightest=abs(rsm["T"] - ref_s),
                             smooth_flux_error=rsm["flux_error"], smooth_nfev=rsm["nfev"]))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/ode1d_validation")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    sharp, smooth, conv = sharp_cases(), smooth_cases(), rtol_convergence()
    write(os.path.join(out, "sharp.csv"), sharp)
    write(os.path.join(out, "smooth.csv"), smooth)
    write(os.path.join(out, "rtol_convergence.csv"), conv)
    mx = lambda rows, k: float(max(r[k] for r in rows))
    summary = dict(
        defaults=dict(method="DOP853", rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL),
        sharp=dict(n=len(sharp), max_err_T_formula=mx(sharp, "err_T_formula"),
                   max_err_t_matching=mx(sharp, "err_t_matching"), max_err_r_matching=mx(sharp, "err_r_matching"),
                   max_flux_error=mx(sharp, "flux_error"), max_range_change_T=mx(sharp, "range_change_T"),
                   max_range_change_t=mx(sharp, "range_change_t")),
        smooth=dict(n=len(smooth), max_rtol_change_T=mx(smooth, "rtol_change_T"),
                    max_range_change_T=mx(smooth, "range_change_T"), max_flux_error=mx(smooth, "flux_error"),
                    intervals={str(r["V_b"]): [r["x_left"], r["x_right"]] for r in smooth}),
        criteria=dict(sharp_T_vs_formula=1e-8, amplitudes_vs_matching=1e-8, flux=1e-9,
                      smooth_rtol_change=1e-9, smooth_range_change=1e-10),
        wall_time_s=time.time() - t0)
    s, m, c = summary["sharp"], summary["smooth"], summary["criteria"]
    summary["pass"] = dict(
        sharp_T_vs_formula=s["max_err_T_formula"] <= c["sharp_T_vs_formula"],
        amplitudes_vs_matching=max(s["max_err_t_matching"], s["max_err_r_matching"]) <= c["amplitudes_vs_matching"],
        flux=max(s["max_flux_error"], m["max_flux_error"]) <= c["flux"],
        smooth_rtol_change=m["max_rtol_change_T"] <= c["smooth_rtol_change"],
        smooth_range_change=m["max_range_change_T"] <= c["smooth_range_change"])
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    write_dir_metadata(out, make_metadata("ode1d", dict(experiment="ode1d_validation",
                                                        energies=EPS_PACKET.tolist()),
                                          wall_time_s=summary["wall_time_s"]))
    print(json.dumps({k: summary[k] for k in ("sharp", "smooth", "pass", "wall_time_s")}, indent=1))
