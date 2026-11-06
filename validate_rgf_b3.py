r"""
Benchmark B3 (docs/PHASE2_METHOD_SELECTION.md §3): one fixed disorder
realization per configuration. RGF (stationary.py) vs an independent direct
sparse Green-function reference (sparse_reference.py), unitarity,
reciprocity, and numerical convergence. Results: docs/RGF_B3_VALIDATION.md.

This is a numerical-method validation of fixed realizations. It says nothing
about the statistics of the disorder-averaged transport.

Predeclared configurations (spec §3 B3; seed 0, L_d = 1):
  A: W = 100, zeta = 0.02    B: W = 100, zeta = 0.05
  C: W = 100, zeta = 0.20    D: W = 200, zeta = 0.05
Five energies per configuration, fixed BEFORE any B3 result existed from the
project's packet-energy convention (spec B4): kx = kx0 + j sigma_k, j = -2..2,
sigma_k = 1/(sqrt(2) sigma_x) = 3.5355, E = (kx^2 + pi^2/L_y^2)/2.
(kx0 = 20, sigma_x = 0.2 -> E = 88.5, 140.5, 204.9, 281.9, 371.4.) The
threshold rule of spec B1/B4 (|E - E_n| >= 1e-3 E) is checked; it moves none.

Base discretization: dx = 1/512, N_y = 63, N_q = 2 N_y + 1 = 127, buffers >= 0.2.
Device: slices x_j on the lattice 1.9 + k dx (the lattice of the slab start),
from 1.9 - n_b dx to 2.9 + n_b dx with the buffer length n_b dx fixed
physically at ceil(buffer * 512)/512 for every dx: 0.2012 (103/512) and
0.4004 (205/512). The disorder region (1.9, 2.9) is physical and identical
for every buffer. Devices of different dx or buffer are nested on one
lattice (the dx = 1/512 slices are a subset of the dx = 1/1024 slices).

Predeclared criteria (unchanged from the spec):
  max|S_RGF - S_sparse| <= 1e-10 (complex, full S)
  U = max(|S^dag S - I|, |S S^dag - I|) <= 1e-10
  Q = max|S - S^T| <= 1e-8 (full), blocks r_L = r_L^T, r_R = r_R^T,
      t_LR = t_RL^T <= 1e-8, probabilities |t_RL[m,n]|^2 = |t_LR[n,m]|^2 <= 1e-8
  |dT_1| <= 2e-4 for N_y 63 -> 127 (and 255 only if zeta = 0.02 fails at 127),
      dx 1/512 -> 1/1024, N_q -> 2 N_q, buffer 0.2 -> 0.4
  same realization: potential samples bitwise equal to a direct
      Potential(cfg, disorder_from_config(cfg, seed)).numpy call
Same-direction transmission symmetry is report-only.

Writes results/rgf_b3_validation/{rgf_sparse_comparison, unitarity,
reciprocity, closed_channel_convergence, dx_convergence,
quadrature_convergence, buffer_convergence, potential_identity}.csv,
summary.json, run_metadata.json.
Usage: python validate_rgf_b3.py [--out DIR] [--cases ABCD]
       python validate_rgf_b3.py --record-tests [--out DIR]   (runs and records the test suites)
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")    # performance only (stationary.py docstring)

import argparse
import csv
import datetime as _dt
import hashlib
import json
import math
import subprocess
import sys
import time

import numpy as np

import smatrix_spec as sm
import sparse_reference as sr
from config import waveguide_config
from disorder import disorder_from_config
from potential import Potential
from provenance import code_fingerprint, environment_info
from stationary import build_device, scattering_matrix

# ----------------------------------------------------------------------
# predeclared protocol (spec §3 B3; do not change)
# ----------------------------------------------------------------------
CASES = {"A": dict(seed=0, W=100.0, Ld=1.0, zeta=0.02), "B": dict(seed=0, W=100.0, Ld=1.0, zeta=0.05),
         "C": dict(seed=0, W=100.0, Ld=1.0, zeta=0.20), "D": dict(seed=0, W=200.0, Ld=1.0, zeta=0.05)}
KX_OFFSETS = (-2, -1, 0, 1, 2)          # in units of sigma_k
THRESHOLD_REL = 1e-3
LIMITS = dict(S_agreement=1e-10, U=1e-10, Q=1e-8, dT1=2e-4)
BASE = dict(dx=1 / 512, Ny=63, Nq=127, buffer=0.2)
VARIANTS = {                                   # one parameter changed at a time (see Ny rows below)
    "dx": dict(BASE, dx=1 / 1024),
    "Nq": dict(BASE, Nq=2 * BASE["Nq"]),
    "buffer": dict(BASE, buffer=0.4),
    "Ny127": dict(BASE, Ny=127, Nq=255),       # criterion: default N_q = 2 N_y + 1 at each N_y
    "Ny127_Nq254": dict(BASE, Ny=127, Nq=254),  # exploratory: pure truncation at fixed N_q = 254 (vs the Nq run)
    "Ny255": dict(BASE, Ny=255, Nq=511),       # conditional (zeta = 0.02 only, if it fails at 127)
}


def passes(value: float, limit: float) -> bool:
    """A criterion is met iff value <= limit (NaN never passes)."""
    return bool(value <= limit)


# ----------------------------------------------------------------------
# configuration, energies, geometry
# ----------------------------------------------------------------------
def case_config(case: str):
    p = CASES[case]
    return waveguide_config(disorder_strength=p["W"], disorder_corr_length=p["zeta"], slab_length=p["Ld"])


def case_potential(case: str):
    """The existing pathway: Potential(cfg, disorder_from_config(cfg, seed))."""
    cfg = case_config(case)
    return cfg, Potential(cfg, disorder_from_config(cfg, seed=CASES[case]["seed"]))


def case_energies(cfg, Ly=1.0):
    from stationary import transverse_basis
    Ep = transverse_basis(Ly, 255).E_perp
    sk = 1.0 / (math.sqrt(2.0) * cfg.sigma_x)
    out = []
    for j in KX_OFFSETS:
        kx = cfg.kx0 + j * sk
        E = 0.5 * (kx ** 2 + (math.pi / Ly) ** 2)
        d = np.abs(E - Ep)
        n = int(np.argmin(d))
        if d[n] < THRESHOLD_REL * E:
            raise RuntimeError(f"energy {E} within {THRESHOLD_REL} E of threshold {n + 1}; rule would apply")
        out.append(dict(j=j, kx=kx, E=E, nearest_threshold_n=n + 1, nearest_threshold_E=float(Ep[n]),
                        rel_dist_to_threshold=float(d[n] / E), n_open_exact=int(np.sum(Ep < E))))
    return out


def buffer_length(buffer):
    """Physical buffer length: the requested length rounded up on the coarse
    (dx = 1/512) lattice, so it is the SAME at every dx (a dx refinement must
    not also change the buffer): 0.2 -> 103/512 = 0.2012, 0.4 -> 205/512 = 0.4004."""
    return math.ceil(buffer / BASE["dx"] - 1e-9) * BASE["dx"]


def geometry(cfg, dx, buffer):
    a, b = cfg.disorder_region
    nb = int(round(buffer_length(buffer) / dx))
    x_start = a - nb * dx
    n_int = int(round((b - a) / dx)) + 2 * nb
    return dict(x_start=x_start, x_end=x_start + n_int * dx, n_slices=n_int + 1, n_buffer_slices=nb,
                buffer_left=nb * dx, buffer_right=nb * dx)


def sha(a) -> str:
    return "sha256:" + hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


class Recorder:
    def __init__(self, f):
        self.f, self.calls = f, []

    def __call__(self, X, Y):
        V = self.f(X, Y)
        self.calls.append(np.array(V, copy=True))
        return V


def direct_reference_samples(case, x_start, dx, n_slices, Ly, Nq):
    """The same physical realization through a FRESH Potential/disorder
    object and independently built coordinates (the 'direct call')."""
    cfg, pot = case_potential(case)
    xs = x_start + dx * np.arange(n_slices)
    y = np.arange(1, Nq + 1) * Ly / (Nq + 1)
    Xg, Yg = np.meshgrid(xs, y, indexing="ij")
    return xs, y, pot.numpy(Xg, Yg)


def compare_arrays(a, b):
    bit = bool(a.shape == b.shape and np.array_equal(a, b))
    d = np.abs(a - b)
    return dict(bitwise_equal=bit, max_abs_diff=float(d.max()),
                max_rel_diff=float(np.max(d / np.maximum(np.abs(b), 1e-300))))


def make_rgf_device(case, cfg, pot, var):
    geo = geometry(cfg, var["dx"], var["buffer"])
    rec = Recorder(pot.numpy)
    dev = build_device(rec, geo["x_start"], geo["x_end"], var["dx"], cfg.Ly, var["Ny"], var["Nq"])
    xs, y, ref = direct_reference_samples(case, geo["x_start"], var["dx"], dev.n_slices, cfg.Ly, var["Nq"])
    cmp = compare_arrays(rec.calls[0], ref)
    cmp.update(coords_bitwise_equal=bool(np.array_equal(xs, dev.xs) and np.array_equal(y, dev.basis.y)))
    return dev, geo, dict(shape=list(ref.shape), potential_hash=sha(ref), x_hash=sha(xs), y_hash=sha(y),
                          V_min=float(ref.min()), V_max=float(ref.max()), V_lead_max=dev.V_lead_max, **cmp), ref


# ----------------------------------------------------------------------
# S-matrix diagnostics
# ----------------------------------------------------------------------
def unitarity(S):
    I = np.eye(S.shape[0])
    a = float(np.max(np.abs(S.conj().T @ S - I)))
    b = float(np.max(np.abs(S @ S.conj().T - I)))
    return dict(U_StS=a, U_SSt=b, U=max(a, b))


def reciprocity(S, n_open):
    p = sm.split_S(S, n_open)
    rL, rR, tLR, tRL = p["r_L"], p["r_R"], p["t_LR"], p["t_RL"]
    P_RL, P_LR = np.abs(tRL) ** 2, np.abs(tLR) ** 2        # P_RL[m, n]: in left n, out right m; P_LR[n, m]: in right m, out left n
    out = dict(Q_full=float(np.max(np.abs(S - S.T))), Q_rL=float(np.max(np.abs(rL - rL.T))),
               Q_rR=float(np.max(np.abs(rR - rR.T))), Q_t=float(np.max(np.abs(tLR - tRL.T))),
               Q_prob=float(np.max(np.abs(P_RL - P_LR.T))),
               same_dir_RL_report_only=float(np.max(np.abs(P_RL - P_RL.T))),
               same_dir_LR_report_only=float(np.max(np.abs(P_LR - P_LR.T))))
    out["Q"] = max(out["Q_full"], out["Q_rL"], out["Q_rR"], out["Q_t"], out["Q_prob"])
    return out


def s_difference(S_a, S_b):
    """Elementwise |S_a - S_b| of the COMPLEX matrices (real and imaginary parts)."""
    S_a, S_b = np.asarray(S_a), np.asarray(S_b)
    if S_a.shape != S_b.shape:
        raise ValueError(f"S shapes differ: {S_a.shape} vs {S_b.shape}")
    return np.abs(S_a - S_b)


def T1(S, n_open):
    t_RL = sm.split_S(S, n_open)["t_RL"]
    return float(np.sum(np.abs(t_RL[:, 0]) ** 2))


# ----------------------------------------------------------------------
# one case
# ----------------------------------------------------------------------
def run_case(case, log=print):
    cfg, pot = case_potential(case)
    p = CASES[case]
    E_list = case_energies(cfg)
    R = dict(identity=[], comparison=[], unitarity=[], reciprocity=[], closed=[], dx=[], quad=[], buffer=[],
             energies=[dict(case=case, **e) for e in E_list], geometry={})
    base = {}                                     # energy index -> base T1 / n_open

    def ident(vname, var, info, geo, extra=None):
        R["identity"].append(dict(case=case, seed=p["seed"], W=p["W"], zeta=p["zeta"], Ld=p["Ld"],
                                  variant=vname, evaluator="RGF device (build_device) vs direct Potential call",
                                  dx=var["dx"], N_y=var["Ny"], N_q=var["Nq"], buffer=var["buffer"],
                                  n_slices=geo["n_slices"], x_start=geo["x_start"], x_end=geo["x_end"],
                                  disorder_region=list(cfg.disorder_region), n_disorder_modes=cfg.disorder_n_modes,
                                  disorder_edge=cfg.disorder_edge, **info, **(extra or {})))

    def run_variant(vname, var):
        dev, geo, info, ref = make_rgf_device(case, cfg, pot, var)
        ident(vname, var, info, geo)
        R["geometry"][vname] = dict(var, **geo, V_lead_max=info["V_lead_max"], potential_hash=info["potential_hash"])
        out = []
        t0 = time.time()
        for i, e in enumerate(E_list):
            r = scattering_matrix(e["E"], dev)
            out.append(dict(r=r, T1=T1(r["S"], r["n_open"]), n_open=r["n_open"], **unitarity(r["S"])))
        log(f"  case {case} variant {vname}: N={dev.n_slices}, N_y={var['Ny']}: {time.time() - t0:.1f}s")
        return dev, geo, ref, out

    # --- base: RGF + sparse reference ---------------------------------
    dev, geo, ref, outs = run_variant("base", BASE)
    for i, (e, o) in enumerate(zip(E_list, outs)):
        r = o["r"]
        t0 = time.time()
        sp_ = sr.sparse_scattering_matrix(e["E"], pot.numpy, geo["x_start"], BASE["dx"], dev.n_slices, cfg.Ly,
                                          BASE["Ny"], BASE["Nq"])
        t_sp = time.time() - t0
        if i == 0:
            R["identity"].append(dict(
                case=case, seed=p["seed"], W=p["W"], zeta=p["zeta"], Ld=p["Ld"], variant="base",
                evaluator="sparse reference (own sampling) vs direct Potential call", dx=BASE["dx"], N_y=BASE["Ny"],
                N_q=BASE["Nq"], buffer=BASE["buffer"], n_slices=dev.n_slices, x_start=geo["x_start"],
                x_end=geo["x_end"], disorder_region=list(cfg.disorder_region), n_disorder_modes=cfg.disorder_n_modes,
                disorder_edge=cfg.disorder_edge, shape=list(ref.shape), potential_hash=sha(sp_["V_samples"]),
                x_hash=sha(sp_["xs"]), y_hash=R["identity"][-1]["y_hash"], V_min=float(ref.min()),
                V_max=float(ref.max()), V_lead_max=dev.V_lead_max,
                coords_bitwise_equal=bool(np.array_equal(sp_["xs"], dev.xs)), **compare_arrays(sp_["V_samples"], ref)))
        S_r, S_s = r["S"], sp_["S"]
        n_open = r["n_open"]
        if sp_["n_open"] != n_open or S_s.shape != S_r.shape:
            raise RuntimeError("RGF and sparse reference disagree on the open-channel count")
        diff = s_difference(S_r, S_s)
        ia = np.unravel_index(int(np.argmax(diff)), diff.shape)
        blocks = {k: float(np.max(np.abs(v - sm.split_S(S_s, n_open)[k]))) for k, v in sm.split_S(S_r, n_open).items()}
        R["comparison"].append(dict(
            case=case, seed=p["seed"], W=p["W"], zeta=p["zeta"], Ld=p["Ld"], kx=e["kx"], E=e["E"], n_open=n_open,
            S_dim=int(S_r.shape[0]), n_slices=dev.n_slices, matrix_order=int(dev.n_slices * BASE["Ny"]),
            max_abs_dS=float(diff.max()), argmax_row=int(ia[0]), argmax_col=int(ia[1]),
            max_abs_dS_rL=blocks["r_L"], max_abs_dS_tLR=blocks["t_LR"], max_abs_dS_tRL=blocks["t_RL"],
            max_abs_dS_rR=blocks["r_R"], T1_RGF=o["T1"], T1_sparse=T1(S_s, n_open), dT1_RGF_sparse=o["T1"] - T1(S_s, n_open),
            sparse_relative_residual=sp_["relative_residual"], sparse_M_nnz=sp_["M_nnz"], sparse_lu_nnz=sp_["lu_nnz"],
            sparse_time_s=t_sp, rgf_min_rcond=r["min_rcond"], limit=LIMITS["S_agreement"],
            passed=passes(float(diff.max()), LIMITS["S_agreement"])))
        for meth, S in (("RGF", S_r), ("sparse", S_s)):
            u, q = unitarity(S), reciprocity(S, n_open)
            R["unitarity"].append(dict(case=case, kx=e["kx"], E=e["E"], method=meth, S_dim=int(S.shape[0]), **u,
                                       limit=LIMITS["U"], passed=passes(u["U"], LIMITS["U"])))
            R["reciprocity"].append(dict(case=case, kx=e["kx"], E=e["E"], method=meth, S_dim=int(S.shape[0]), **q,
                                         limit=LIMITS["Q"], passed=passes(q["Q"], LIMITS["Q"])))
        base[i] = dict(T1=o["T1"], n_open=n_open, U=o["U"])

    # --- refinements: one parameter at a time vs the base T1 -----------
    def rows_for(vname, var, role, ref_name="base", ref_vals=None, ref_var=None, target=None):
        _, _, _, vo = run_variant(vname, var)
        rows = []
        for i, (e, o) in enumerate(zip(E_list, vo)):
            b = ref_vals[i] if ref_vals is not None else base[i]
            d = o["T1"] - b["T1"]
            rows.append(dict(
                case=case, seed=p["seed"], W=p["W"], zeta=p["zeta"], Ld=p["Ld"], kx=e["kx"], E=e["E"], role=role,
                reference=ref_name, variant=vname, T1_reference=b["T1"], T1_variant=o["T1"], delta_T1=d,
                abs_delta_T1=abs(d), limit=LIMITS["dT1"], passed=passes(abs(d), LIMITS["dT1"]),
                n_open_reference=b["n_open"], n_open_variant=o["n_open"], U_variant=o["U"]))
        return rows, vo

    # dx
    rows, _ = rows_for("dx", VARIANTS["dx"], "criterion")
    for r in rows:
        r.update(dx_reference=BASE["dx"], dx_variant=VARIANTS["dx"]["dx"])
    R["dx"] += rows
    # Nq
    rows, vq = rows_for("Nq", VARIANTS["Nq"], "criterion")
    for r in rows:
        r.update(Nq_reference=BASE["Nq"], Nq_variant=VARIANTS["Nq"]["Nq"])
    R["quad"] += rows
    nq_vals = [dict(T1=o["T1"], n_open=o["n_open"]) for o in vq]
    # buffer
    geo_b = geometry(cfg, BASE["dx"], 0.4)
    rows, _ = rows_for("buffer", VARIANTS["buffer"], "criterion")
    for r in rows:
        r.update(buffer_reference=geometry(cfg, BASE["dx"], 0.2)["buffer_left"], buffer_variant=geo_b["buffer_left"])
    R["buffer"] += rows
    # closed channels
    def closed(rows, Ny_ref, Ny_var):
        for r in rows:
            r.update(Ny_reference=Ny_ref, Ny_variant=Ny_var, closed_reference=Ny_ref - r["n_open_reference"],
                     closed_variant=Ny_var - r["n_open_variant"])
        return rows
    rows, v127 = rows_for("Ny127", VARIANTS["Ny127"], "criterion")
    R["closed"] += closed(rows, BASE["Ny"], 127)
    v127_vals = [dict(T1=o["T1"], n_open=o["n_open"]) for o in v127]
    rows, _ = rows_for("Ny127_Nq254", VARIANTS["Ny127_Nq254"], "exploratory_fixed_Nq254", ref_name="Nq(63,254)",
                       ref_vals=nq_vals)
    R["closed"] += closed(rows, BASE["Ny"], 127)
    if p["zeta"] == 0.02 and any(not r["passed"] for r in R["closed"] if r["role"] == "criterion"):
        log(f"  case {case}: |dT1| criterion failed at N_y = 127 -> running N_y = 255 (spec B3 rule)")
        rows, _ = rows_for("Ny255", VARIANTS["Ny255"], "conditional_extension_127_to_255", ref_name="Ny127",
                           ref_vals=v127_vals)
        R["closed"] += closed(rows, 127, 255)
        rows, _ = rows_for("Ny255", VARIANTS["Ny255"], "conditional_extension_63_to_255", ref_name="base")
        R["closed"] += closed(rows, BASE["Ny"], 255)
    return R


# ----------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------
def _crit(rows, key, limit, extra=("kx", "E")):
    vals = np.array([r[key] for r in rows], float)
    i = int(np.nanargmax(vals)) if np.all(np.isfinite(vals)) else int(np.argmax(~np.isfinite(vals)))
    cases = sorted({r["case"] for r in rows})
    per = {c: float(max(r[key] for r in rows if r["case"] == c)) for c in cases}
    return dict(limit=limit, value=float(vals[i]), argmax={"case": rows[i]["case"], **{k: rows[i][k] for k in extra},
                                                         **({"method": rows[i]["method"]} if "method" in rows[i] else {})},
                n=len(rows), per_case_max=per, per_case_pass={c: passes(v, limit) for c, v in per.items()},
                passed=bool(all(passes(v, limit) for v in per.values())))


def summarize(R):
    s = {}
    s["S_agreement_RGF_vs_sparse"] = _crit(R["comparison"], "max_abs_dS", LIMITS["S_agreement"])
    s["S_agreement_RGF_vs_sparse"]["S_dims"] = sorted({r["S_dim"] for r in R["comparison"]})
    for meth in ("RGF", "sparse"):
        rows = [r for r in R["unitarity"] if r["method"] == meth]
        s[f"unitarity_{meth}"] = dict(_crit(rows, "U", LIMITS["U"]), max_StS=float(max(r["U_StS"] for r in rows)),
                                      max_SSt=float(max(r["U_SSt"] for r in rows)))
        rows = [r for r in R["reciprocity"] if r["method"] == meth]
        s[f"reciprocity_{meth}"] = dict(
            full=_crit(rows, "Q_full", LIMITS["Q"]), rL=_crit(rows, "Q_rL", LIMITS["Q"]),
            rR=_crit(rows, "Q_rR", LIMITS["Q"]), tLR_vs_tRL_T=_crit(rows, "Q_t", LIMITS["Q"]),
            probability=_crit(rows, "Q_prob", LIMITS["Q"]),
            report_only_same_direction=dict(
                max_RL=float(max(r["same_dir_RL_report_only"] for r in rows)),
                max_LR=float(max(r["same_dir_LR_report_only"] for r in rows)),
                note="not a criterion: same-direction symmetry does not follow from time reversal (spec §2.7.8)"))
        s[f"reciprocity_{meth}"]["passed"] = all(s[f"reciprocity_{meth}"][k]["passed"]
                                                 for k in ("full", "rL", "rR", "tLR_vs_tRL_T", "probability"))
    crit = lambda rows: [r for r in rows if r["role"] == "criterion"]
    for name, rows in (("closed_channel_Ny63_vs_Ny127", crit(R["closed"])), ("dx_1_512_vs_1_1024", crit(R["dx"])),
                       ("Nq_127_vs_254", crit(R["quad"])), ("buffer_0.2_vs_0.4", crit(R["buffer"]))):
        s[name] = _crit(rows, "abs_delta_T1", LIMITS["dT1"])
        s[name]["max_abs_delta_T1_per_case_energy"] = {f"{r['case']}:{r['kx']:.4f}": r["abs_delta_T1"] for r in rows}
        s[name]["rows_failed"] = [dict(case=r["case"], kx=r["kx"], E=r["E"], abs_delta_T1=r["abs_delta_T1"])
                                  for r in rows if not r["passed"]]
    cl = s["closed_channel_Ny63_vs_Ny127"]
    cl["closed_channels_retained_range_Ny63"] = [int(min(r["closed_reference"] for r in crit(R["closed"]))),
                                                 int(max(r["closed_reference"] for r in crit(R["closed"])))]
    cl["closed_channels_retained_range_Ny127"] = [int(min(r["closed_variant"] for r in crit(R["closed"]))),
                                                  int(max(r["closed_variant"] for r in crit(R["closed"])))]
    ext = [r for r in R["closed"] if r["role"].startswith("conditional")]
    s["closed_channel_conditional_Ny255"] = dict(
        triggered=bool(ext), rule="N_y = 255 is run only if the zeta = 0.02 case fails at N_y = 127",
        rows=[dict(case=r["case"], kx=r["kx"], E=r["E"], role=r["role"], delta_T1=r["delta_T1"],
                   abs_delta_T1=r["abs_delta_T1"], passed=r["passed"]) for r in ext],
        max_abs_delta_127_to_255=float(max([r["abs_delta_T1"] for r in ext if r["role"].endswith("127_to_255")],
                                           default=float("nan"))),
        max_abs_delta_63_to_255=float(max([r["abs_delta_T1"] for r in ext if r["role"].endswith("63_to_255")],
                                          default=float("nan"))),
        status="exploratory-labelled extension of the protocol; the criterion result above is not replaced")
    expl = [r for r in R["closed"] if r["role"].startswith("exploratory")]
    s["exploratory_fixed_Nq254_Ny63_vs_127"] = dict(
        max_abs_delta_T1=float(max(r["abs_delta_T1"] for r in expl)),
        note="exploratory, not a criterion: isolates the transverse truncation (N_q = 254 at both N_y)")
    ident = R["identity"]
    s["potential_identity"] = dict(
        passed=bool(all(r["bitwise_equal"] and r["coords_bitwise_equal"] for r in ident)), n_checks=len(ident),
        max_abs_diff=float(max(r["max_abs_diff"] for r in ident)),
        max_V_lead_on_first_lead_slices=float(max(r["V_lead_max"] for r in ident)),
        note="bitwise equality of the sampled potential with a direct Potential(cfg, disorder_from_config(cfg, seed)).numpy call")
    keys = [k for k, v in s.items() if isinstance(v, dict) and "passed" in v and k != "closed_channel_conditional_Ny255"]
    s["criteria_passed"] = {k: bool(s[k]["passed"]) for k in keys}
    s["all_B3_criteria_passed"] = bool(all(s["criteria_passed"].values()))
    s["overall_B3_status"] = "PASS" if s["all_B3_criteria_passed"] else "FAIL (see criteria_passed)"
    return s


# ----------------------------------------------------------------------
def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)


def _pytest(args):
    p = subprocess.run([sys.executable, "-m", "pytest", "-q"] + args, capture_output=True, text=True)
    last = [l for l in p.stdout.strip().splitlines() if l.strip()][-1]
    return dict(command="python -m pytest -q " + " ".join(args), exit_status=p.returncode, summary=last)


def record_tests(out):
    path = os.path.join(out, "run_metadata.json")
    with open(path) as f:
        m = json.load(f)
    m["post_run_tests"] = dict(
        recorded_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        runs=[_pytest(["tests/test_rgf_b3.py"]), _pytest(["tests"])],
        note="executed after the artifacts were generated; written by validate_rgf_b3.py --record-tests")
    with open(path, "w") as f:
        json.dump(m, f, indent=2)
    print(json.dumps(m["post_run_tests"], indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rgf_b3_validation")
    ap.add_argument("--cases", default="ABCD")
    ap.add_argument("--record-tests", action="store_true")
    a = ap.parse_args()
    if a.record_tests:
        record_tests(a.out)
        sys.exit(0)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    R = {}
    for case in a.cases:
        print(f"case {case}: {CASES[case]}", flush=True)
        r = run_case(case, log=lambda m: print(m, flush=True))
        for k, v in r.items():
            if k in ("energies",) or isinstance(v, list):
                R.setdefault(k, []).extend(v)
            else:
                R.setdefault(k, {})[case] = v
    names = dict(comparison="rgf_sparse_comparison", unitarity="unitarity", reciprocity="reciprocity",
                 closed="closed_channel_convergence", dx="dx_convergence", quad="quadrature_convergence",
                 buffer="buffer_convergence", identity="potential_identity")
    for k, fn in names.items():
        write_csv(os.path.join(a.out, fn + ".csv"), R[k])
    summary = summarize(R)
    summary.update(wall_time_s=time.time() - t0, limits=LIMITS,
                   scope="numerical-method validation of fixed disorder realizations; not a statistical result",
                   cases=CASES, base_discretization=BASE)
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    env = environment_info()
    import scipy
    meta = dict(
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        command=" ".join([os.path.basename(sys.executable)] + sys.argv),
        source_commit=env.get("git_commit"), working_tree_dirty=env.get("git_dirty"),
        python=env.get("python"), numpy=np.__version__, scipy=scipy.__version__, platform=env.get("platform"),
        machine=env.get("machine"), cpu_count=env.get("cpu_count"),
        openblas_num_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
        rgf_implementation="stationary.scattering_matrix -> stationary.rgf_endpoint_blocks (left-to-right block recursion)",
        sparse_reference="sparse_reference.sparse_scattering_matrix: scipy.sparse.linalg.splu (SuperLU, COLAMD), "
                         "direct solves for 2 N_y boundary right-hand sides; no recursion; no import of stationary.py",
        shared_conventions="smatrix_spec.lead_modes / fisher_lee / assemble_S; same discretized Hamiltonian",
        cases={c: dict(CASES[c], disorder_n_modes=case_config(c).disorder_n_modes, disorder_edge=case_config(c).disorder_edge,
                       disorder_region=list(case_config(c).disorder_region),
                       disorder_params_hash=sha(np.concatenate([disorder_from_config(case_config(c), seed=CASES[c]["seed"]).qx,
                                                                disorder_from_config(case_config(c), seed=CASES[c]["seed"]).qy,
                                                                disorder_from_config(case_config(c), seed=CASES[c]["seed"]).phi])))
               for c in a.cases},
        energies=R["energies"], base_discretization=BASE, variants=VARIANTS, geometry=R["geometry"],
        energy_rule="kx = kx0 + j sigma_k (j = -2..2), sigma_k = 1/(sqrt(2) sigma_x), E = (kx^2 + pi^2/L_y^2)/2; "
                    "fixed before any B3 result; threshold rule |E - E_n| >= 1e-3 E checked (moves none)",
        limits=LIMITS, wall_time_s=summary["wall_time_s"], **code_fingerprint("rgf_b3"))
    with open(os.path.join(a.out, "run_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps({k: summary[k] for k in ("criteria_passed", "overall_B3_status", "wall_time_s")}, indent=1))
