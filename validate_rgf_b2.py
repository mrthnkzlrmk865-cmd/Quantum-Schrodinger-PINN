r"""
Benchmark B2 (docs/PHASE2_METHOD_SELECTION.md §3): stationary RGF
(stationary.py) versus the independent 1D ODE reference (ode1d.py) for the
smooth rectangular barrier of waveguide_config, uniform in y.
Results: docs/RGF_B2_VALIDATION.md.

Same potential in both calculations: Potential(cfg).barrier_numpy, with
cfg = waveguide_config(barrier_height=V_b), V_b in {150, 250}, width 0.1,
edge 0.01, centre 2.4, no disorder. The RGF samples it at its slice
positions (bitwise check below), the ODE evaluates it at its integration
points (ode1d.potential_barrier_callable). Nothing else is shared.

Energy convention: the ODE solves -1/2 phi'' + V_b phi = eps phi; the RGF is
run at the total energy E = eps + E_perp_1, E_perp_1 = pi^2/(2 L_y^2) (the
same double the RGF basis uses). Probabilities are compared:
T_RGF(E) = sum_m |(t_RL)_{m1}|^2 vs T_ODE(eps).

Energies: kx = linspace(4.09, 35.91, 60), eps = kx^2/2. A total energy
closer than 1e-3 E to a transverse threshold (spec B1/B4 rule) is moved to
E_n (1 +/- 2e-3) on its own side; the same adjusted eps is used by the ODE.
(For the 60 B2 energies no point is that close; the closest is 1.9e-3 E.)

Binding criteria (spec §3 B2, unchanged):
  |T_RGF - T_ODE| <= 2e-4 at dx = 1/512, all 60 energies, both V_b;
  off-diagonal T^{L->R}_{m!=1 <- 1} <= 1e-12; U <= 1e-10; Q = max|S - S^T| <= 1e-8;
  ODE reference: flux <= 1e-9, rtol 1e-10 -> 1e-12 change <= 1e-9,
  clean range doubled <= 1e-10, sharp barrier vs closed form <= 1e-8.
Measured, not assumed: convergence order from dx = 1/256, 1/512, 1/1024.
Provisional (reported): Richardson |dT| <= 2e-5, applied only where the
measured order is in [1.8, 2.2].

Separate part (--td): TD (split-step, production settings, committed runs
results/transport/barrier/Vb*_W0.npz, re-run live for reproducibility) vs
the ODE averaged over the packet spectrum P(kx) (spec B4 protocol).

Writes results/rgf_b2_validation/{rgf_ode_comparison.csv, ode_reference.csv,
summary.json, run_metadata.json, td_energy_average.json, *.png}.
Usage: python validate_rgf_b2.py [--out DIR] [--no-td] [--no-plots]
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")    # performance only (stationary.py docstring)

import argparse
import csv
import datetime as _dt
import json
import sys
import time

import numpy as np

import smatrix_spec as sm
from analytic import rectangular_barrier_transmission
from config import waveguide_config
from ode1d import (DEFAULT_ATOL, DEFAULT_RTOL, clean_interval, potential_barrier_callable,
                   rectangular_barrier, solve_scattering)
from potential import Potential
from provenance import code_fingerprint, environment_info
from stationary import build_device, scattering_matrix, transverse_basis

# ----------------------------------------------------------------------
# predeclared parameters and thresholds (spec §3 B2 / B4; do not change)
# ----------------------------------------------------------------------
VB_LIST = (150.0, 250.0)
DXS = (1 / 256, 1 / 512, 1 / 1024)
BINDING_DX = 1 / 512
LY = 1.0
NY = 63
X_START, X_END = 1.5, 3.5               # device (buffers around the barrier at 2.35..2.45)
X_START_WIDE, X_END_WIDE = 1.0, 4.0     # buffer check at the binding dx
KX = np.linspace(4.09, 35.91, 60)
THRESHOLD_REL = 1e-3                    # |E - E_n| >= 1e-3 E (spec B1/B4)
THRESHOLD_SHIFT = 2e-3                  # deterministic move: E_n (1 +/- 2e-3)
ODE_RANGE_TOL = 1e-12                   # clean_interval tolerance (validate_ode1d.py)
ODE_LOOSE = 1e-10

LIMITS = dict(dT=2e-4, mixing=1e-12, U=1e-10, Q=1e-8)
PROVISIONAL = dict(richardson_dT=2e-5, order_window=(1.8, 2.2))
ODE_LIMITS = dict(sharp_T_vs_formula=1e-8, flux=1e-9, smooth_rtol_change=1e-9, smooth_range_change=1e-10)
TD_LIMITS = dict(interval=1e-3, P_in=1e-4, quad_halving=2e-4)
ORDER_FLOOR = 1e-8                      # |error| below this (10x ODE accuracy) makes an order unreliable


# ----------------------------------------------------------------------
# energies
# ----------------------------------------------------------------------
def e_perp():
    """Transverse threshold energies, the same doubles the RGF uses."""
    return transverse_basis(LY, NY).E_perp


def total_energy(eps, Ep1):
    return eps + Ep1


def adjust_energy(E, E_perp, rel=THRESHOLD_REL, shift=THRESHOLD_SHIFT):
    """Deterministic threshold rule. If |E - E_n| < rel*E for the nearest
    threshold E_n, return E_n (1 + s shift), s = +1 if E >= E_n else -1
    (the point stays on its own side, so the number of open channels does not
    change). Otherwise E is returned unchanged (bitwise)."""
    E = float(E)
    d = np.abs(E - E_perp)
    n = int(np.argmin(d))
    out = dict(E_original=E, nearest_threshold_n=n + 1, nearest_threshold_E=float(E_perp[n]),
               dist_original=float(d[n]), rel_dist_original=float(d[n] / E), adjusted=False, E=E)
    if d[n] < rel * E:
        s = 1.0 if E >= E_perp[n] else -1.0
        Ea = float(E_perp[n] * (1.0 + s * shift))
        if np.any(np.abs(Ea - E_perp) < rel * Ea):
            raise RuntimeError(f"threshold rule failed at E = {E}")
        out.update(adjusted=True, E=Ea)
    dn = np.abs(out["E"] - E_perp)
    out.update(dist=float(dn.min()), rel_dist=float(dn.min() / out["E"]),
               n_open=int(np.sum(E_perp < out["E"])))
    return out


def energy_grid(Ep=None):
    Ep = e_perp() if Ep is None else Ep
    Ep1 = float(Ep[0])
    rows = []
    for i, kx in enumerate(KX):
        eps0 = 0.5 * float(kx) ** 2
        a = adjust_energy(total_energy(eps0, Ep1), Ep)
        eps = eps0 if not a["adjusted"] else a["E"] - Ep1
        rows.append(dict(i=i, kx=float(kx), eps_original=eps0, eps=eps, E_total=total_energy(eps, Ep1), **{k: v for k, v in a.items() if k != "E"}))
    return rows


# ----------------------------------------------------------------------
# RGF side
# ----------------------------------------------------------------------
class RecordingPotential:
    """Wraps Potential.barrier_numpy as the V(X, Y) callable of build_device
    and records every array it returns (to check what the RGF sampled)."""

    def __init__(self, pot):
        self.pot = pot
        self.calls = []

    def __call__(self, X, Y):
        V = self.pot.barrier_numpy(X)
        self.calls.append((np.array(X, copy=True), np.array(V, copy=True)))
        return V


def rgf_device(pot, dx, x_start=X_START, x_end=X_END):
    rec = RecordingPotential(pot)
    dev = build_device(rec, x_start, x_end, dx, LY, NY)
    Xs, Vs = rec.calls[0]                              # device slices x quadrature points
    ref = pot.barrier_numpy(dev.xs)
    bitwise = bool(np.array_equal(Xs[:, 0], dev.xs) and np.all(Xs == Xs[:, :1])
                   and all(np.array_equal(Vs[:, l], ref) for l in range(Vs.shape[1])))
    V = dev.V
    off = V - np.einsum("jnn->jn", V)[:, :, None] * np.eye(NY)[None]
    chk = dict(n_slices=dev.n_slices, n_quad=dev.basis.n_quad, samples_bitwise_equal_barrier_numpy=bitwise,
               max_offdiag_V_block=float(np.max(np.abs(off))),
               max_diag_V_block_minus_Vx=float(np.max(np.abs(np.einsum("jnn->jn", V) - ref[:, None]))),
               max_V=float(ref.max()), V_lead_max=dev.V_lead_max,
               V_at_device_ends=[float(ref[0]), float(ref[-1])])
    return dev, chk


def rgf_point(E, dev):
    r = scattering_matrix(E, dev)
    t_RL, t_LR, r_L = r["t_RL"], r["t_LR"], r["r_L"]
    S = r["S"]
    I = np.eye(S.shape[0])
    PT, PTb = np.abs(t_RL) ** 2, np.abs(t_LR) ** 2
    off = lambda P: P - np.diag(np.diag(P))
    T = float(np.sum(PT[:, 0]))
    R = float(np.sum(np.abs(r_L[:, 0]) ** 2))
    return dict(
        n_open=r["n_open"], T_RGF=T, R_RGF=R, T_RGF_11=float(PT[0, 0]),
        T_RGF_R_to_L_ch1=float(np.sum(PTb[:, 0])),
        unitarity_residual=float(max(np.max(np.abs(S.conj().T @ S - I)), np.max(np.abs(S @ S.conj().T - I)))),
        reciprocity_residual=float(sm.reciprocity_residual(S)),
        mode_reciprocity_residual=float(sm.mode_reciprocity_residual(t_RL, t_LR)),
        channel_mixing=float(max(np.max(off(PT)[:, 0]), np.max(off(PTb)[:, 0]))),
        channel_mixing_all_channels=float(max(np.max(off(PT)), np.max(off(PTb)))),
        reflection_mixing_ch1=float(np.max(off(np.abs(r_L) ** 2)[:, 0])),
        flux_residual_ch1=abs(T + R - 1.0), min_rcond=r["min_rcond"])


# ----------------------------------------------------------------------
# ODE side
# ----------------------------------------------------------------------
def ode_setup(cfg):
    V = potential_barrier_callable(cfg)
    x1 = cfg.barrier_center - cfg.barrier_width / 2
    x2 = cfg.barrier_center + cfg.barrier_width / 2
    a, b = clean_interval(V, x1, x2, tol=ODE_RANGE_TOL)
    return V, x1, x2, a, b


def ode_point(eps, cfg, setup):
    V, x1, x2, a, b = setup
    base = solve_scattering(eps, V, a, b)
    loose = solve_scattering(eps, V, a, b, rtol=ODE_LOOSE, atol=ODE_LOOSE)
    wide = solve_scattering(eps, V, a - (x1 - a), b + (b - x2))
    Vs, bp = rectangular_barrier(cfg.barrier_height, 0.0, cfg.barrier_width)
    sharp = solve_scattering(eps, Vs, 0.0, cfg.barrier_width, bp)
    Tf = rectangular_barrier_transmission(np.array([eps]), cfg.barrier_height, cfg.barrier_width)[0]
    return dict(T_ODE=base["T"], R_ODE=base["R"], ODE_flux_error=base["flux_error"], ODE_nfev=base["nfev"],
                ODE_x_left=a, ODE_x_right=b, ODE_rtol=base["rtol"], ODE_atol=base["atol"],
                ODE_method=base["method"], ODE_rtol_change_T=abs(loose["T"] - base["T"]),
                ODE_range_change_T=abs(wide["T"] - base["T"]),
                ODE_sharp_err_vs_formula=abs(sharp["T"] - Tf), ODE_sharp_flux_error=sharp["flux_error"])


# ----------------------------------------------------------------------
# convergence order
# ----------------------------------------------------------------------
def observed_order(e_coarse, e_fine, ratio=2.0, floor=ORDER_FLOOR):
    """log(|e_coarse|/|e_fine|)/log(ratio); NaN if either error is below the
    floor (not resolvable against the reference accuracy)."""
    if abs(e_coarse) < floor or abs(e_fine) < floor:
        return float("nan")
    return float(np.log(abs(e_coarse) / abs(e_fine)) / np.log(ratio))


def richardson(T_coarse, T_fine, p=2.0, ratio=2.0):
    return T_fine + (T_fine - T_coarse) / (ratio ** p - 1.0)


def in_window(p, w=PROVISIONAL["order_window"]):
    return bool(np.isfinite(p) and w[0] <= p <= w[1])


def regime(eps, Vb):
    return "below" if eps < 0.9 * Vb else ("near" if eps <= 1.1 * Vb else "above")


# ----------------------------------------------------------------------
# main RGF vs ODE benchmark
# ----------------------------------------------------------------------
def run_benchmark(grid, log=print):
    rows, ode_rows, dev_checks, buffer_rows = [], [], {}, []
    for Vb in VB_LIST:
        cfg = waveguide_config(barrier_height=Vb)
        pot = Potential(cfg)
        setup = ode_setup(cfg)
        t0 = time.time()
        ode = {}
        for g in grid:
            ode[g["i"]] = ode_point(g["eps"], cfg, setup)
            ode_rows.append(dict(Vb=Vb, kx=g["kx"], eps=g["eps"], **ode[g["i"]]))
        log(f"  ODE V_b={Vb}: {len(grid)} energies in {time.time() - t0:.1f}s")
        rgf = {}
        for dx in DXS:
            dev, chk = rgf_device(pot, dx)
            dev_checks[f"Vb={Vb},dx=1/{round(1 / dx)}"] = chk
            t0 = time.time()
            for g in grid:
                rgf[(g["i"], dx)] = rgf_point(g["E_total"], dev)
            log(f"  RGF V_b={Vb} dx=1/{round(1 / dx)} N={dev.n_slices}: {time.time() - t0:.1f}s")
        devw, chkw = rgf_device(pot, BINDING_DX, X_START_WIDE, X_END_WIDE)
        dev_checks[f"Vb={Vb},dx=1/{round(1 / BINDING_DX)},wide"] = chkw
        for g in grid:
            Tw = rgf_point(g["E_total"], devw)["T_RGF"]
            buffer_rows.append(dict(Vb=Vb, kx=g["kx"], buffer_change_T=abs(Tw - rgf[(g["i"], BINDING_DX)]["T_RGF"])))
        for g in grid:
            o = ode[g["i"]]
            e = {dx: rgf[(g["i"], dx)]["T_RGF"] - o["T_ODE"] for dx in DXS}
            p_c = observed_order(e[DXS[0]], e[DXS[1]])
            p_f = observed_order(e[DXS[1]], e[DXS[2]])
            T3 = [rgf[(g["i"], dx)]["T_RGF"] for dx in DXS]
            p_self = observed_order(T3[0] - T3[1], T3[1] - T3[2], floor=0.0)
            rich_ok = in_window(p_c) and in_window(p_f)
            Tr = richardson(T3[1], T3[2]) if rich_ok else float("nan")
            for dx, p in zip(DXS, (float("nan"), p_c, p_f)):
                r = rgf[(g["i"], dx)]
                rows.append(dict(
                    Vb=Vb, width=cfg.barrier_width, edge=cfg.barrier_edge, kx=g["kx"], epsilon=g["eps"],
                    E_total=g["E_total"], energy_adjusted=g["adjusted"], regime=regime(g["eps"], Vb), dx=dx,
                    N_y=NY, N_q=dev_checks[f"Vb={Vb},dx=1/{round(1 / dx)}"]["n_quad"],
                    number_of_slices=dev_checks[f"Vb={Vb},dx=1/{round(1 / dx)}"]["n_slices"],
                    number_of_open_channels=r["n_open"], T_RGF=r["T_RGF"], R_RGF=r["R_RGF"],
                    T_ODE=o["T_ODE"], R_ODE=o["R_ODE"], signed_T_error=e[dx], absolute_T_error=abs(e[dx]),
                    relative_T_error=abs(e[dx]) / o["T_ODE"] if o["T_ODE"] >= 1e-6 else float("nan"),
                    observed_convergence_order=p, order_coarse_pair=p_c, order_fine_pair=p_f,
                    order_self_three_grid=p_self, richardson_applied=rich_ok,
                    T_richardson=Tr, richardson_abs_error=abs(Tr - o["T_ODE"]) if rich_ok else float("nan"),
                    unitarity_residual=r["unitarity_residual"], reciprocity_residual=r["reciprocity_residual"],
                    mode_reciprocity_residual=r["mode_reciprocity_residual"],
                    channel_mixing=r["channel_mixing"], channel_mixing_all_channels=r["channel_mixing_all_channels"],
                    reflection_mixing_ch1=r["reflection_mixing_ch1"], T_RGF_11=r["T_RGF_11"],
                    T_RGF_R_to_L_ch1=r["T_RGF_R_to_L_ch1"], flux_residual_ch1=r["flux_residual_ch1"],
                    min_rcond=r["min_rcond"], ODE_flux_error=o["ODE_flux_error"], ODE_nfev=o["ODE_nfev"]))
    return rows, ode_rows, dev_checks, buffer_rows


def _stats(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if x.size == 0:
        return dict(n=0)
    q = np.percentile(x, [0, 25, 50, 75, 100])
    return dict(n=int(x.size), min=float(q[0]), q25=float(q[1]), median=float(q[2]), q75=float(q[3]),
                max=float(q[4]))


def summarize(rows, ode_rows, dev_checks, buffer_rows, grid):
    s = dict(per_barrier={})
    for Vb in VB_LIST:
        R = [r for r in rows if r["Vb"] == Vb]
        per_dx = {}
        for dx in DXS:
            rr = [r for r in R if r["dx"] == dx]
            err = np.array([r["absolute_T_error"] for r in rr])
            i = int(np.argmax(err))
            per_dx[f"1/{round(1 / dx)}"] = dict(
                N=rr[0]["number_of_slices"], n_energies=len(rr), max_abs_T_error=float(err.max()),
                rms_abs_T_error=float(np.sqrt(np.mean(err ** 2))),
                argmax=dict(kx=rr[i]["kx"], epsilon=rr[i]["epsilon"], E_total=rr[i]["E_total"]),
                max_by_regime={g: float(max([r["absolute_T_error"] for r in rr if r["regime"] == g], default=np.nan))
                               for g in ("below", "near", "above")},
                max_unitarity=float(max(r["unitarity_residual"] for r in rr)),
                max_reciprocity=float(max(r["reciprocity_residual"] for r in rr)),
                max_mode_reciprocity=float(max(r["mode_reciprocity_residual"] for r in rr)),
                max_channel_mixing=float(max(r["channel_mixing"] for r in rr)),
                max_channel_mixing_all_channels=float(max(r["channel_mixing_all_channels"] for r in rr)),
                max_flux_residual_ch1=float(max(r["flux_residual_ch1"] for r in rr)),
                min_rcond=float(min(r["min_rcond"] for r in rr)))
        e512 = [r for r in R if r["dx"] == DXS[1]]
        e1024 = [r for r in R if r["dx"] == DXS[2]]
        mx = [per_dx[f"1/{round(1 / dx)}"]["max_abs_T_error"] for dx in DXS]
        rich = [r["richardson_abs_error"] for r in e1024 if r["richardson_applied"]]
        O = [o for o in ode_rows if o["Vb"] == Vb]
        bf = [b["buffer_change_T"] for b in buffer_rows if b["Vb"] == Vb]
        s["per_barrier"][str(Vb)] = dict(
            per_dx=per_dx,
            binding_dx_max_abs_T_error=per_dx["1/512"]["max_abs_T_error"],
            binding_dx_argmax=per_dx["1/512"]["argmax"],
            order_coarse_pair_256_512=_stats([r["order_coarse_pair"] for r in e512]),
            order_fine_pair_512_1024=_stats([r["order_fine_pair"] for r in e1024]),
            order_self_three_grid=_stats([r["order_self_three_grid"] for r in e1024]),
            n_order_nan_coarse=int(sum(not np.isfinite(r["order_coarse_pair"]) for r in e512)),
            n_order_nan_fine=int(sum(not np.isfinite(r["order_fine_pair"]) for r in e1024)),
            n_sign_change_512_1024=int(sum(np.sign(a["signed_T_error"]) != np.sign(b["signed_T_error"])
                                           for a, b in zip(e512, e1024))),
            order_of_max_errors=[float(np.log2(mx[0] / mx[1])), float(np.log2(mx[1] / mx[2]))],
            richardson=dict(n_applied=len(rich), n_energies=len(e1024),
                            max_abs_error=float(max(rich)) if rich else None,
                            provisional_pass_where_applied=bool(rich and max(rich) <= PROVISIONAL["richardson_dT"])),
            ode=dict(max_flux_error=float(max(o["ODE_flux_error"] for o in O)),
                     max_rtol_change_T=float(max(o["ODE_rtol_change_T"] for o in O)),
                     max_range_change_T=float(max(o["ODE_range_change_T"] for o in O)),
                     max_sharp_err_vs_formula=float(max(o["ODE_sharp_err_vs_formula"] for o in O)),
                     max_sharp_flux_error=float(max(o["ODE_sharp_flux_error"] for o in O)),
                     nfev_range=[int(min(o["ODE_nfev"] for o in O)), int(max(o["ODE_nfev"] for o in O))],
                     interval=[O[0]["ODE_x_left"], O[0]["ODE_x_right"]],
                     method=O[0]["ODE_method"], rtol=O[0]["ODE_rtol"], atol=O[0]["ODE_atol"]),
            rgf_buffer_check=dict(device_wide=[X_START_WIDE, X_END_WIDE], dx=BINDING_DX,
                                  max_change_T=float(max(bf))))
    PB = s["per_barrier"]
    allr = rows
    crit = dict(
        A_binding_dT_dx512=dict(limit=LIMITS["dT"], value={k: v["binding_dx_max_abs_T_error"] for k, v in PB.items()}),
        C_channel_mixing=dict(limit=LIMITS["mixing"], value=float(max(r["channel_mixing"] for r in allr))),
        D_unitarity=dict(limit=LIMITS["U"], value=float(max(r["unitarity_residual"] for r in allr))),
        E_reciprocity=dict(limit=LIMITS["Q"], value=float(max(r["reciprocity_residual"] for r in allr))),
        F_ode_flux=dict(limit=ODE_LIMITS["flux"], value=float(max(max(v["ode"]["max_flux_error"],
                                                                      v["ode"]["max_sharp_flux_error"])
                                                                  for v in PB.values()))),
        F_ode_rtol_change=dict(limit=ODE_LIMITS["smooth_rtol_change"],
                               value=float(max(v["ode"]["max_rtol_change_T"] for v in PB.values()))),
        F_ode_range_change=dict(limit=ODE_LIMITS["smooth_range_change"],
                                value=float(max(v["ode"]["max_range_change_T"] for v in PB.values()))),
        F_ode_sharp_vs_formula=dict(limit=ODE_LIMITS["sharp_T_vs_formula"],
                                    value=float(max(v["ode"]["max_sharp_err_vs_formula"] for v in PB.values()))))
    for k, c in crit.items():
        v = c["value"]
        c["pass"] = bool(max(v.values()) <= c["limit"]) if isinstance(v, dict) else bool(v <= c["limit"])
    crit["B_convergence_order"] = dict(
        status="measured and reported (not a pass/fail criterion)",
        richardson="applied only where both per-energy orders (vs ODE) are in [1.8, 2.2]; "
                   "provisional |dT| <= 2e-5 reported for those energies only")
    sb = all(dc["samples_bitwise_equal_barrier_numpy"] for dc in dev_checks.values())
    s.update(
        criteria=crit, limits=LIMITS, provisional=PROVISIONAL, ode_limits=ODE_LIMITS,
        potential_checks=dict(all_samples_bitwise_equal=sb, per_device=dev_checks),
        energy_grid=dict(n=len(grid), n_adjusted=int(sum(g["adjusted"] for g in grid)),
                         min_rel_dist_to_threshold=float(min(g["rel_dist"] for g in grid)),
                         rule=f"|E - E_n| < {THRESHOLD_REL} E -> E_n (1 +/- {THRESHOLD_SHIFT}), same side",
                         n_open_range=[min(g["n_open"] for g in grid), max(g["n_open"] for g in grid)]),
        rgf_vs_ode_all_binding_pass=bool(all(c["pass"] for c in crit.values() if "pass" in c) and sb))
    return s


# ----------------------------------------------------------------------
# TD vs energy-averaged ODE (separate part; spec B2 last bullet with B4 protocol)
# ----------------------------------------------------------------------
TD_DIR = "results/transport/barrier"
TD_T_END, TD_N_FRAMES = 1.5, 41
DKX = 0.025
N_SIGMA = 4.5


def packet_weights(kx, kx0, sigma_x):
    return sigma_x / np.sqrt(np.pi) * np.exp(-((kx - kx0) ** 2) * sigma_x ** 2)


def kx_window(kx0, sigma_x, dkx=DKX, n_sigma=N_SIGMA):
    sk = 1.0 / (np.sqrt(2.0) * sigma_x)
    lo, hi = kx0 - n_sigma * sk, kx0 + n_sigma * sk
    n = int(round((hi - lo) / dkx)) + 1
    return np.linspace(lo, hi, n)


def energy_average(kx, vals, kx0, sigma_x):
    """Trapezoid integral of P(kx) vals(kx) on a uniform kx grid; also returns int P."""
    P = packet_weights(kx, kx0, sigma_x)
    return float(np.trapezoid(P * vals, kx)), float(np.trapezoid(P, kx))


def interval_distance(x, lo, hi):
    return 0.0 if lo <= x <= hi else float(min(abs(x - lo), abs(x - hi)))


def run_td(log=print):
    from provenance import fingerprint, read_npz_metadata, solver_config_dict
    from transport import run_transport
    from experiment_disorder import transport_spec
    out = {}
    for Vb in VB_LIST:
        name = f"Vb{int(Vb)}_W0"
        path = os.path.join(TD_DIR, f"{name}.npz")
        cfg = waveguide_config(disorder_strength=0.0, barrier_height=Vb)
        d = np.load(path, allow_pickle=False)
        meta = read_npz_metadata(path)          # None: these files predate the provenance metadata
        fp_ok = None if meta is None else bool(
            meta.get("fingerprint") == fingerprint(transport_spec(cfg, 1, TD_T_END, TD_N_FRAMES)))
        stored = json.loads(str(d["config_json"]))
        solver_fields = json.loads(json.dumps(solver_config_dict(cfg)))   # PINN-only fields excluded
        cfg_ok = bool(all(stored.get(k) == v for k, v in solver_fields.items()))
        T_td, R_td, P_td = float(d["T"][-1, 0]), float(d["R"][-1, 0]), float(d["P_in"][-1, 0])
        t0 = time.time()
        live = run_transport(cfg, [0], t_end=TD_T_END, n_frames=TD_N_FRAMES)
        rep = dict(max_abs_diff_T_R_P=float(max(np.max(np.abs(live[k] - d[k])) for k in ("T", "R", "P_in"))),
                   bitwise=bool(all(np.array_equal(live[k], d[k]) for k in ("T", "R", "P_in"))),
                   wall_time_s=time.time() - t0)
        # informative TD refinements (not the production result, nothing written to results/transport)
        refine = {}
        for lab, over in (("dx_1/128", dict(nx=int(round(cfg.Lx * 128)) - 1)),
                          ("dx_1/256", dict(nx=int(round(cfg.Lx * 256)) - 1)), ("dt_1e-4", dict(dt=1e-4))):
            rr = run_transport(cfg.replace(**over), [0], t_end=TD_T_END, n_frames=TD_N_FRAMES)
            refine[lab] = dict(T=float(rr["T"][-1, 0]), R=float(rr["R"][-1, 0]), P_in=float(rr["P_in"][-1, 0]),
                               dT_vs_production=float(rr["T"][-1, 0] - T_td), wall_time_s=rr["wall_time_s"])
        # ODE energy average over the packet spectrum
        setup = ode_setup(cfg)
        kf = kx_window(cfg.kx0, cfg.sigma_x, DKX / 2)
        kc = kx_window(cfg.kx0, cfg.sigma_x, DKX)
        assert np.allclose(kf[::2], kc, rtol=0, atol=1e-12)
        E_tot = 0.5 * (kf ** 2 + (np.pi / cfg.Ly) ** 2)            # E(kx) of spec B4
        eps = E_tot - 0.5 * (np.pi / cfg.Ly) ** 2                   # = kx^2/2 (ODE energy)
        t0 = time.time()
        sol = [solve_scattering(e, setup[0], setup[3], setup[4]) for e in eps]
        Tf = np.array([r["T"] for r in sol])
        Rf = np.array([r["R"] for r in sol])
        avT_f, normP_f = energy_average(kf, Tf, cfg.kx0, cfg.sigma_x)
        avR_f, _ = energy_average(kf, Rf, cfg.kx0, cfg.sigma_x)
        avT_c, normP_c = energy_average(kf[::2], Tf[::2], cfg.kx0, cfg.sigma_x)
        avR_c, _ = energy_average(kf[::2], Rf[::2], cfg.kx0, cfg.sigma_x)
        T0 = solve_scattering(0.5 * cfg.kx0 ** 2, setup[0], setup[3], setup[4])["T"]
        dT = interval_distance(avT_c, T_td, T_td + P_td)
        dR = interval_distance(avR_c, R_td, R_td + P_td)
        out[str(Vb)] = dict(
            td_file=path, td_npz_has_provenance_metadata=meta is not None,
            td_fingerprint_matches_production_spec=fp_ok, td_config_equals_waveguide_config=cfg_ok,
            td_times_end=float(d["times"][-1]), td_n_frames=int(d["times"].size),
            td_settings=dict(nx=cfg.nx, ny=cfg.ny, dt=cfg.dt, Lx=cfg.Lx, t_end=TD_T_END, kx0=cfg.kx0,
                             sigma_x=cfg.sigma_x, x0=cfg.x0, cap_width=cfg.cap_width, cap_strength=cfg.cap_strength,
                             T_region="x > 2.9 + absorbed right", R_region="x < 1.9 + absorbed left"),
            T_TD=T_td, R_TD=R_td, P_in_TD=P_td, td_live_rerun=rep, td_refinements_informative=refine,
            kx_window=[float(kc[0]), float(kc[-1])], n_kx=int(kc.size), dkx=float(kc[1] - kc[0]),
            int_P=normP_c, int_P_half_step=normP_f,
            window_truncation_bound=1.0 - normP_c,     # |change of <T>, <R>| from the finite window <= 1 - int P
            avg_T_ODE=avT_c, avg_R_ODE=avR_c, avg_T_ODE_half_step=avT_f, avg_R_ODE_half_step=avR_f,
            quad_halving_change_T=abs(avT_f - avT_c), quad_halving_change_R=abs(avR_f - avR_c),
            T_ODE_central_energy_not_used=T0,
            ode_max_flux_error=float(max(r["flux_error"] for r in sol)),
            dist_T=dT, dist_R=dR,
            pass_P_in=bool(P_td <= TD_LIMITS["P_in"]),
            pass_quad=bool(max(abs(avT_f - avT_c), abs(avR_f - avR_c)) <= TD_LIMITS["quad_halving"]),
            pass_T=bool(dT <= TD_LIMITS["interval"]), pass_R=bool(dR <= TD_LIMITS["interval"]),
            ode_wall_time_s=time.time() - t0)
        v = out[str(Vb)]
        v["pass"] = bool(cfg_ok and rep["bitwise"] and v["pass_P_in"] and v["pass_quad"] and v["pass_T"] and v["pass_R"])
        log(f"  TD V_b={Vb}: T_TD={T_td:.6f} <T>_ODE={avT_c:.6f} dist={dT:.2e}; R dist={dR:.2e}; "
            f"rerun bitwise={rep['bitwise']}")
    return dict(limits=TD_LIMITS, protocol="spec §3 B4: P(kx) = sigma_x/sqrt(pi) exp[-(kx-kx0)^2 sigma_x^2], "
                "E(kx) = (kx^2 + pi^2/L_y^2)/2, uniform kx over kx0 +/- 4.5/(sqrt(2) sigma_x), dkx = 0.025, "
                "trapezoid; halving check dkx/2; TD interval [T_TD, T_TD + P_in]",
                barriers=out, all_pass=bool(all(v["pass"] for v in out.values())))


# ----------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------
COLORS = {DXS[0]: "#2a78d6", DXS[1]: "#eb6834", DXS[2]: "#1baf7a"}
MARKERS = {DXS[0]: "o", DXS[1]: "s", DXS[2]: "^"}


def make_plots(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted = "#0b0b0b", "#52514e"
    plt.rcParams.update({"axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})
    # 1: |dT| vs epsilon for the three dx
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, Vb in zip(axs, VB_LIST):
        for dx in DXS:
            rr = [r for r in rows if r["Vb"] == Vb and r["dx"] == dx]
            ax.semilogy([r["epsilon"] for r in rr], [max(r["absolute_T_error"], 1e-16) for r in rr],
                        marker=MARKERS[dx], ms=4, lw=1.2, color=COLORS[dx], label=f"dx = 1/{round(1 / dx)}")
        ax.axhline(LIMITS["dT"], color=ink, lw=1, ls="--", label="2e-4 (binding at dx = 1/512)")
        ax.axvline(Vb, color=muted, lw=0.8, ls=":")
        ax.set_title(f"V_b = {Vb:g}", color=ink)
        ax.set_xlabel("longitudinal energy ε = kx²/2")
    axs[0].set_ylabel("|T_RGF − T_ODE|")
    axs[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "error_vs_energy.png"), dpi=150)
    plt.close(fig)
    # 2: |dT| vs dx, one line per energy, with slope-2 guide
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, Vb in zip(axs, VB_LIST):
        eps_vals = sorted({r["epsilon"] for r in rows if r["Vb"] == Vb})
        cmap = plt.get_cmap("Blues")
        for k, e in enumerate(eps_vals):
            rr = sorted([r for r in rows if r["Vb"] == Vb and r["epsilon"] == e], key=lambda r: r["dx"])
            ax.loglog([r["dx"] for r in rr], [max(r["absolute_T_error"], 1e-16) for r in rr], lw=0.8,
                      color=cmap(0.35 + 0.6 * k / len(eps_vals)), marker=".", ms=3)
        dd = np.array(DXS)
        ax.loglog(dd, 2e-4 * (dd / DXS[1]) ** 2, color=ink, lw=1.2, ls="--", label="slope 2 (guide)")
        ax.axhline(LIMITS["dT"], color="#eb6834", lw=1, ls=":", label="2e-4")
        ax.set_title(f"V_b = {Vb:g} (60 energies; darker = higher ε)", color=ink)
        ax.set_xlabel("dx")
    axs[0].set_ylabel("|T_RGF − T_ODE|")
    axs[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "error_vs_dx.png"), dpi=150)
    plt.close(fig)
    # 3: observed order vs energy
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, Vb in zip(axs, VB_LIST):
        r512 = [r for r in rows if r["Vb"] == Vb and r["dx"] == DXS[1]]
        r1024 = [r for r in rows if r["Vb"] == Vb and r["dx"] == DXS[2]]
        ax.plot([r["epsilon"] for r in r512], [r["order_coarse_pair"] for r in r512], "o", ms=4,
                color=COLORS[DXS[0]], label="1/256 → 1/512")
        ax.plot([r["epsilon"] for r in r1024], [r["order_fine_pair"] for r in r1024], "^", ms=4,
                color=COLORS[DXS[2]], label="1/512 → 1/1024")
        ax.axhspan(1.8, 2.2, color="#e4e3df", zorder=0, label="[1.8, 2.2]")
        ax.axvline(Vb, color=muted, lw=0.8, ls=":")
        ax.set_title(f"V_b = {Vb:g}", color=ink)
        ax.set_xlabel("ε")
    axs[0].set_ylabel("observed order (vs ODE)")
    axs[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "order_vs_energy.png"), dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rgf_b2_validation")
    ap.add_argument("--no-td", action="store_true", help="skip the separate TD energy-averaged comparison")
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    Ep = e_perp()
    grid = energy_grid(Ep)
    write_csv(os.path.join(a.out, "energy_grid.csv"), grid)
    rows, ode_rows, dev_checks, buffer_rows = run_benchmark(grid, log=lambda m: print(m, flush=True))
    write_csv(os.path.join(a.out, "rgf_ode_comparison.csv"), rows)
    write_csv(os.path.join(a.out, "ode_reference.csv"), ode_rows)
    write_csv(os.path.join(a.out, "rgf_buffer_check.csv"), buffer_rows)
    summary = summarize(rows, ode_rows, dev_checks, buffer_rows, grid)
    td = None
    if not a.no_td:
        td = run_td(log=lambda m: print(m, flush=True))
        with open(os.path.join(a.out, "td_energy_average.json"), "w") as f:
            json.dump(td, f, indent=2)
    summary["td_energy_averaged_comparison"] = (
        "not run" if td is None else dict(all_pass=td["all_pass"], file="td_energy_average.json",
                                          per_barrier={k: dict(dist_T=v["dist_T"], dist_R=v["dist_R"],
                                                               pass_all=v["pass"]) for k, v in td["barriers"].items()}))
    summary["overall_B2_status"] = (
        "PASS" if summary["rgf_vs_ode_all_binding_pass"] and td is not None and td["all_pass"] else
        "FAIL" if not summary["rgf_vs_ode_all_binding_pass"] or (td is not None and not td["all_pass"]) else
        "INCOMPLETE (TD comparison not run)")
    summary["conventions"] = dict(
        energy="E_total = epsilon + E_perp_1, E_perp_1 = pi^2/(2 L_y^2) (stationary.transverse_basis value)",
        rgf_reference_planes=f"convention P: first and last device slices x = {X_START} and {X_END} (spec §2.7.5)",
        ode_reference_planes="incident/reflected at x_left, transmitted at x_right of the clean_interval (ode1d)",
        compared="probabilities T, R (independent of the reference planes)")
    summary["wall_time_s"] = time.time() - t0
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    if not a.no_plots:
        make_plots(rows, a.out)
    env = environment_info()
    import scipy
    meta = dict(
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        command=" ".join([os.path.basename(sys.executable)] + sys.argv),
        source_commit=env.get("git_commit"), working_tree_dirty=env.get("git_dirty"),
        python=env.get("python"), numpy=np.__version__, scipy=scipy.__version__,
        platform=env.get("platform"), machine=env.get("machine"), cpu_count=env.get("cpu_count"),
        openblas_num_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
        parameters=dict(Vb=list(VB_LIST), width=0.1, edge=0.01, barrier_center=2.4, Ly=LY, N_y=NY,
                        N_q=2 * NY + 1, dxs=list(DXS), binding_dx=BINDING_DX, device=[X_START, X_END],
                        device_wide_check=[X_START_WIDE, X_END_WIDE], kx=KX.tolist(),
                        threshold_rule=dict(rel=THRESHOLD_REL, shift=THRESHOLD_SHIFT),
                        ode=dict(method="DOP853", rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL, loose=ODE_LOOSE,
                                 clean_interval_tol=ODE_RANGE_TOL),
                        order_floor=ORDER_FLOOR, td=dict(t_end=TD_T_END, n_frames=TD_N_FRAMES, dkx=DKX,
                                                         n_sigma=N_SIGMA)),
        limits=dict(binding=LIMITS, provisional=PROVISIONAL, ode=ODE_LIMITS, td=TD_LIMITS),
        wall_time_s=summary["wall_time_s"], **code_fingerprint("rgf_b2"))
    with open(os.path.join(a.out, "run_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps({k: summary[k] for k in ("criteria", "overall_B2_status", "wall_time_s")}, indent=1))
