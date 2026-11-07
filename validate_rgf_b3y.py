r"""
B3Y driver: independent fine-grid validation of the RGF at the predeclared failed-B3X points, with the Tier 3 exact-lattice reference.
Protocol (fixed before any result): docs/RGF_B3Y_VALIDATION.md Part I and docs/rgf_b3y_protocol.json (revision 2). Results: Part II.

Methods (each keeps its own implementation; nothing is shared between them except the physical potential function):
  M1  stationary.scattering_matrix                    RGF block recursion + smatrix_spec Fisher-Lee
  M2  sparse_reference.sparse_scattering_matrix       SuperLU (COLAMD) Green blocks + smatrix_spec Fisher-Lee
  M3  sparse_wavefunction.wavefunction_scattering_matrix   own lead data, cosine-sum blocks, SuperLU (MMD_AT_PLUS_A) + 1 refinement, S from the state
  M4  sparse_exact_lattice                            unmerged operator, extended-precision residuals, cancellation-free lead data

Tier 3 error decomposition (M2 and M3 ONLY; M1 has no explicit matrix and receives none):
  S_double - S_exact = (S_hybrid - S_exact) + (S_fl - S_hybrid) + (S_double - S_fl)
  S_fl      the method's OWN rounded matrix, double lead data in matrix, sources and extraction, solved exactly in extended precision
  S_hybrid  the same rounded bulk operator with ONLY the two lead (corner) terms replaced by the exact extended-precision ones;
            exact-lead sources and extraction

Usage:  python validate_rgf_b3y.py [--out DIR] [--cache DIR] [--stages a,b,...]
        python validate_rgf_b3y.py --synthetic --out DIR --cache DIR          (no B3Y energy point is used)
        python validate_rgf_b3y.py --record-tests --out DIR
The driver refuses to write into the directories of earlier validation results.
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
import pickle
import platform
import subprocess
import sys
import time
import types

import numpy as np

import sparse_exact_lattice as xp
import sparse_reference as sr
import sparse_wavefunction as sw
import validate_rgf_b3 as b3
import validate_rgf_b3x as b3x
from provenance import code_fingerprint, environment_info
from smatrix_spec import fisher_lee as _fisher_lee_unused  # noqa: F401  (imported so the dependency is explicit)
from stationary import build_device, scattering_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
PROTOCOL_PATH = os.path.join(HERE, "docs", "rgf_b3y_protocol.json")
B3X_CSV = os.path.join(HERE, "results", "rgf_b3x_validation", "per_energy_comparison.csv")
with open(PROTOCOL_PATH) as _f:
    PROTOCOL = json.load(_f)
LIMITS = PROTOCOL["limits"]
XPG = PROTOCOL["xp_gates"]
DISC = PROTOCOL["discretization"]
FORBIDDEN_OUTPUT_DIRS = tuple(os.path.normpath(os.path.join(HERE, "results", d)) for d in (
    "rgf_b1_validation", "rgf_b2_validation", "rgf_b3_validation", "rgf_b3x_validation", "ode1d_validation", "validation", "reproduction"))
EPS_MACH = 2.220446049250313e-16
LD, CLD = xp.LD, xp.CLD
Z = types.SimpleNamespace   # tiny record type


def passes(v, lim):
    return bool(v <= lim)              # NaN never passes


# ======================================================================
# problems (the physical device family); the synthetic problem is a small stand-in used only for tests and the synthetic driver run
# ======================================================================
class RealProblem:
    """The predeclared B3Y problem: B3X configurations, seed 0, N_y = 63, N_q = 127."""
    name = "B3Y"
    Ly, Ny, Nq = 1.0, DISC["N_y"], DISC["N_q"]
    inv_base = 512
    invs = tuple(round(1 / d) for d in DISC["dx_levels"])           # (512, 1024, 2048)
    buffer_req, buffer_ctrl_req = 0.2, 0.4
    seed = PROTOCOL["realization"]["seed"]

    def __init__(self):
        self._cfg = {}
        points = PROTOCOL["points"]
        self.points = {p["id"]: dict(p) for p in points}
        S = PROTOCOL["point_sets"]
        self.sets = {k: list(S[k]) for k in ("tier1", "tier2_dx2048", "tier3_xp_1024_2048", "tier3_xp_512_validation", "buffer_control")}
        self.sets_map = dict(tier1=self.sets["tier1"], tier2=self.sets["tier2_dx2048"], tier3_fine=self.sets["tier3_xp_1024_2048"],
                             tier3_validation=self.sets["tier3_xp_512_validation"], buffer=self.sets["buffer_control"])

    def potential(self, cid):
        if cid not in self._cfg:
            self._cfg[cid] = b3x.config_potential(cid)
        return self._cfg[cid]

    def region(self, cid):
        return tuple(self.potential(cid)[0].disorder_region)

    def direct(self, cid, x_start, dx, n):
        return b3x.direct_samples(cid, x_start, dx, n, self.Ly, self.Nq)

    def geometry(self, cid, inv, buffer):
        cfg, _ = self.potential(cid)
        return b3.geometry(cfg, 1.0 / inv, buffer)


def _syn_numpy(X, Y):
    """Smooth, y-dependent synthetic potential that mixes channels; clean (|V| < 1e-9) at the lead slices."""
    return (20.0 * np.exp(-((X - 1.5) / 0.1) ** 2) * (1.0 + 0.5 * np.cos(2 * np.pi * Y))
            + 5.0 * np.exp(-((X - 1.5) / 0.15) ** 2) * np.sin(3 * np.pi * Y))


class SyntheticProblem:
    """Small stand-in (no B3Y energy, no B3Y realization): N_y = 7, N_q = 15, dx = 1/32, 1/64, 1/128."""
    name = "SYNTHETIC"
    Ly, Ny, Nq = 1.0, 7, 15
    inv_base = 32
    invs = (32, 64, 128)
    buffer_req, buffer_ctrl_req = 0.25, 0.5
    seed = None

    def __init__(self):
        E_a, E_b = 60.0, 30.0
        Ep = 0.5 * (np.arange(1, self.Ny + 1) * np.pi / self.Ly) ** 2
        self.points = {}
        for pid, E in (("SYN_A@E60", E_a), ("SYN_B@E30", E_b)):
            self.points[pid] = dict(id=pid, config="SYN", E=E, n_open=int(np.sum(Ep < E)), kx=float("nan"), role="synthetic")
        a, b = "SYN_A@E60", "SYN_B@E30"
        self.sets_map = dict(tier1=[a, b], tier2=[a], tier3_fine=[a], tier3_validation=[a, b], buffer=[a])
        self.sets = self.sets_map

    def potential(self, cid):
        cfg = types.SimpleNamespace(Ly=self.Ly, disorder_region=(1.0, 2.0))
        return cfg, types.SimpleNamespace(numpy=lambda X, Y: np.asarray(_syn_numpy(X, Y), dtype=float))

    def region(self, cid):
        return (1.0, 2.0)

    def direct(self, cid, x_start, dx, n):
        xs = x_start + dx * np.arange(n)
        y = np.arange(1, self.Nq + 1) * self.Ly / (self.Nq + 1)
        Xg, Yg = np.meshgrid(xs, y, indexing="ij")
        return _syn_numpy(Xg, Yg)

    def geometry(self, cid, inv, buffer):
        a, b = self.region(cid)
        dx0 = 1.0 / self.inv_base
        blen = math.ceil(buffer / dx0 - 1e-9) * dx0
        dx = 1.0 / inv
        nb = int(round(blen / dx))
        x_start = a - nb * dx
        n_int = int(round((b - a) / dx)) + 2 * nb
        return dict(x_start=x_start, x_end=x_start + n_int * dx, n_slices=n_int + 1, n_buffer_slices=nb, buffer_left=nb * dx, buffer_right=nb * dx)


# ======================================================================
# method-specific building blocks used by the decomposition (mirrors of the methods' OWN source / extraction conventions)
# ======================================================================
def wf_source_mirror(t, lam, a, o):
    """sparse_wavefunction: s = _source_sign * (-t) * (1/lambda - lambda) * a with _source_sign = 1.0 (same operation order, bitwise)."""
    return 1.0 * (-t) * (1.0 / lam[o] - lam[o]) * a


def wf_extract(c1, cN, o, sg, a):
    """The wave-function extraction of sparse_wavefunction.wavefunction_S_from_blocks, operation for operation; works in any precision.
    c1, cN: (N_y, 2 n_open) boundary slices of the solution columns [left source | right source]; sg = sqrt(Gamma_open)[:, None]; a = Gamma^{-1/2}."""
    no = o.size
    r_L = sg * (c1[o][:, :no] - np.diag(a))
    t_RL = sg * cN[o][:, :no]
    t_LR = sg * c1[o][:, no:]
    r_R = sg * (cN[o][:, no:] - np.diag(a))
    return np.block([[r_L, t_LR], [t_RL, r_R]])


def fl_extract(G11o, GN1o, G1No, GNNo, sg, one_j):
    """smatrix_spec.fisher_lee + assemble_S on open-channel blocks (same operation order; one_j = 1j or CLD(1j)); sg = sqrt(Gamma_open)."""
    no = sg.shape[0]
    f = lambda G: one_j * sg[:, None] * G * sg[None, :]
    I = np.eye(no)
    return np.block([[f(G11o) - I, f(G1No)], [f(GN1o), f(GNNo) - I]])


def _blocks(X, N, ny):
    return X.reshape(N, ny, X.shape[1])


def _term_dict(dec):
    return dict(bulk_rounding=dec["bulk_rounding"], lead_data=dec["lead_data"], solver=dec["solver"])


def _finish_decomposition(S_double, S_exact, S_fl, S_hyb, ex_ind, fl_ind, hy_ind, fl_gate, hy_gate, extra):
    dec = xp.error_decomposition(S_double, S_exact, S_fl, S_hyb)
    e_comb = xp.combined_error_indicator(ex_ind, fl_ind, hy_ind)
    dom = xp.representation_dominance(_term_dict(dec), e_comb)
    return dict(terms=dec, e_ind_exact=ex_ind, e_ind_fl=fl_ind, e_ind_hybrid=hy_ind, e_ind_combined=e_comb, dominance=dom,
                refinement_fl=fl_gate, refinement_hybrid=hy_gate, S_fl=S_fl, S_hybrid=S_hyb, **extra)


def decompose_m3(E, wf, dx, lead, S_exact, ex_ind):
    """Tier 3 decomposition of the wave-function method (M3) from ITS OWN rounded matrix, double lead data, source and extraction."""
    V, Ep = wf["V_blocks"], wf["E_perp"]
    N, ny = V.shape[0], V.shape[1]
    t = sw.hopping(dx)
    lam, gam, op = wf["lam"], wf["gamma"], wf["open"]
    o = np.where(op)[0]
    no = o.size
    M = sw.assemble(E, V, Ep, dx, lam)                                              # the method's rounded matrix
    a = wf["a"]
    src_mirror = wf_source_mirror(t, lam, a, o)
    src_bitwise = bool(np.array_equal(src_mirror, wf["source"]))
    B = np.zeros((N * ny, 2 * no), dtype=complex)
    B[o, np.arange(no)] = wf["source"]
    B[(N - 1) * ny + o, no + np.arange(no)] = wf["source"]
    # extraction mirror reproduces the method's own S from ITS c1, cN (bitwise)
    S_mirror = wf_extract(wf["c1"], wf["cN"], o, np.sqrt(gam[o])[:, None], a)
    ext_bitwise = bool(np.array_equal(S_mirror, wf["S"]))
    sg_d, a_d = np.sqrt(gam[o]).astype(CLD)[:, None], a.astype(CLD)
    ex_d = lambda X: wf_extract(_blocks(X, N, ny)[0], _blocks(X, N, ny)[-1], o, sg_d, a_d)
    fl = xp.solve_matrix_extended(M, B.astype(CLD), permc_spec="MMD_AT_PLUS_A", extract=ex_d)
    S_fl = ex_d(fl["X"])
    if not np.array_equal(lead["open"], op):
        raise ValueError("open-channel masks of the method and of the exact lattice differ: no decomposition")
    sg_x = np.sqrt(lead["gamma"][o]).astype(CLD)[:, None]
    a_x = (lead["gamma"][o] ** LD(-0.5)).astype(CLD)
    ex_x = lambda X: wf_extract(_blocks(X, N, ny)[0], _blocks(X, N, ny)[-1], o, sg_x, a_x)
    Bx = xp.build_sources(lead, N, ny).reshape(N * ny, 2 * no)
    hy = xp.solve_matrix_extended(M, Bx, permc_spec="MMD_AT_PLUS_A", extract=ex_x,
                                  corner_swap=dict(t=float(t), ny=ny, N=N, lam_double=lam, lam_ext=lead["lam"]))
    S_hyb = ex_x(hy["X"])
    return _finish_decomposition(wf["S"], S_exact, S_fl, S_hyb, ex_ind, fl["error_indicator"], hy["error_indicator"],
                                 fl["refinement_gate"], hy["refinement_gate"],
                                 dict(method="M3", mirror_source_bitwise=src_bitwise, mirror_extraction_bitwise=ext_bitwise,
                                      kappa1_fl=fl["kappa1"], kappa1_hybrid=hy["kappa1"]))


def decompose_m2(E, gb, E_perp, dx, lead, S_exact, ex_ind):
    """Tier 3 decomposition of the Green-block method (M2): its OWN matrix (sr.assemble_M), double lead data, unit sources, smatrix_spec Fisher-Lee."""
    V = gb["V_blocks"]
    N, ny = V.shape[0], V.shape[1]
    modes = gb["modes"]
    op, gam = modes["open"], modes["gamma"]
    o = np.where(op)[0]
    no = o.size
    t = xp.hopping(dx)
    M = sr.assemble_M(E, V, E_perp, dx, modes["sigma"])
    G = gb["G"]
    ix = np.ix_(o, o)
    S_mirror = fl_extract(G["G11"][ix], G["GN1"][ix], G["G1N"][ix], G["GNN"][ix], np.sqrt(gam[o]), 1j)
    ext_bitwise = bool(np.array_equal(S_mirror, gb["S"]))
    B = np.zeros((N, ny, 2 * no), dtype=complex)                                    # unit sources on the open channels only (all that the extraction uses)
    B[0, o, np.arange(no)] = 1.0
    B[N - 1, o, no + np.arange(no)] = 1.0
    sg_d = np.sqrt(gam[o]).astype(CLD)

    def make(extract_sg):
        def ex(X):
            Xb = _blocks(X, N, ny)
            return fl_extract(Xb[0][o][:, :no], Xb[-1][o][:, :no], Xb[0][o][:, no:], Xb[-1][o][:, no:], extract_sg, CLD(1j))
        return ex
    ex_d = make(sg_d)
    fl = xp.solve_matrix_extended(M, B.reshape(N * ny, 2 * no).astype(CLD), permc_spec="COLAMD", extract=ex_d)
    S_fl = ex_d(fl["X"])
    if not np.array_equal(lead["open"], op):
        raise ValueError("open-channel masks of the method and of the exact lattice differ: no decomposition")
    ex_x = make(np.sqrt(lead["gamma"][o]).astype(CLD))
    hy = xp.solve_matrix_extended(M, xp.build_unit_sources(lead, N, ny).reshape(N * ny, 2 * no), permc_spec="COLAMD", extract=ex_x,
                                  corner_swap=dict(t=float(t), ny=ny, N=N, lam_double=modes["lam"], lam_ext=lead["lam"]))
    S_hyb = ex_x(hy["X"])
    return _finish_decomposition(gb["S"], S_exact, S_fl, S_hyb, ex_ind, fl["error_indicator"], hy["error_indicator"],
                                 fl["refinement_gate"], hy["refinement_gate"],
                                 dict(method="M2", mirror_source_bitwise=True, mirror_extraction_bitwise=ext_bitwise,
                                      kappa1_fl=fl["kappa1"], kappa1_hybrid=hy["kappa1"]))


def decomposition_row(base, d):
    t = d["terms"]
    return dict(base, method=d["method"], total=t["total"], bulk_rounding=t["bulk_rounding"], lead_data=t["lead_data"], solver=t["solver"],
                representation_total=t["representation_total"], sum_of_term_norms=t.get("sum_of_term_norms"),
                identity_residual_three_term=t["identity_residual_three_term"], identity_residual_two_term=t["identity_residual_two_term"],
                e_ind_exact=d["e_ind_exact"], e_ind_fl=d["e_ind_fl"], e_ind_hybrid=d["e_ind_hybrid"], e_ind_combined=d["e_ind_combined"],
                kappa1_fl=d["kappa1_fl"], kappa1_hybrid=d["kappa1_hybrid"],
                gate_fl_passed=d["refinement_fl"]["passed"], gate_hybrid_passed=d["refinement_hybrid"]["passed"],
                resolved=d["dominance"]["resolved"], dominant=d["dominance"]["dominant"], dominance_supported=d["dominance"]["supported"],
                dominates_each_term_individually=d["dominance"].get("dominates_each_term_individually"),
                bulk_fraction_of_term_budget=d["dominance"].get("bulk_fraction_of_term_budget"), dominance_reason=d["dominance"]["reason"],
                mirror_source_bitwise=d["mirror_source_bitwise"], mirror_extraction_bitwise=d["mirror_extraction_bitwise"])


# ======================================================================
# one case: M1, M2, M3
# ======================================================================
class DeviceCache:
    def __init__(self, prob):
        self.prob, self.d = prob, {}

    def get(self, cid, inv, buffer):
        key = (cid, inv, buffer)
        if key not in self.d:
            if len(self.d) >= 2:
                self.d.pop(next(iter(self.d)))
            cfg, pot = self.prob.potential(cid)
            geo = self.prob.geometry(cid, inv, buffer)
            dx = 1.0 / inv
            rec = b3.Recorder(pot.numpy)
            try:
                dev = build_device(rec, geo["x_start"], geo["x_end"], dx, self.prob.Ly, self.prob.Ny, self.prob.Nq)
                err = None
            except ValueError as e:                                                  # lead-tail guard (P4)
                dev, err = None, str(e)
            ref = self.prob.direct(cid, geo["x_start"], dx, geo["n_slices"] if dev is None else dev.n_slices)
            cmp = b3.compare_arrays(rec.calls[0], ref) if dev is not None else dict(bitwise_equal=False, max_abs_diff=float("nan"), max_rel_diff=float("nan"))
            xs = geo["x_start"] + dx * np.arange(ref.shape[0])
            info = dict(shape=list(ref.shape), potential_hash=b3.sha(ref), x_hash=b3.sha(xs), V_min=float(ref.min()), V_max=float(ref.max()),
                        V_lead_max=float(max(np.max(np.abs(ref[0])), np.max(np.abs(ref[-1])))) if dev is None else dev.V_lead_max,
                        lead_guard_error=err, **cmp)
            self.d[key] = (dev, geo, ref, info, pot)
        return self.d[key]


def _T1s(S, n):
    return b3.T1(S, n)


def rho_diag(dS, E, Ep, dx, N, n_open):
    """Diagnostic only (never classifies): max|dS| / [(N-1) eps / sin(k dx)] for the lowest-k open channel."""
    x = (E - Ep[n_open - 1]) * dx * dx
    sk = math.sqrt(max(x * (2.0 - x), 1e-300))
    return dS / ((N - 1) * EPS_MACH / sk)


def run_methods(prob, devs, pt, inv, buffer, keep_S=True):
    """M1, M2, M3 for one point and resolution. Returns a flat row (+ the three complex S matrices)."""
    cid, E = pt["config"], pt["E"]
    dev, geo, ref, info, pot = devs.get(cid, inv, buffer)
    dx = 1.0 / inv
    row = dict(point=pt["id"], config=cid, E=E, kx=pt.get("kx", float("nan")), inv_dx=inv, dx=dx, buffer_requested=buffer,
               buffer_length=geo["buffer_left"], n_buffer_slices=geo["n_buffer_slices"], x_start=geo["x_start"], x_end=geo["x_end"],
               N_y=prob.Ny, N_q=prob.Nq, V_lead_max=info["V_lead_max"], potential_hash=info["potential_hash"], x_hash=info["x_hash"],
               P1_rgf_device_bitwise=info["bitwise_equal"], P4_lead_guard_ok=bool(info["lead_guard_error"] is None and info["V_lead_max"] <= LIMITS["lead_tail_guard"]))
    if dev is None:
        row.update(status="NOT EVALUABLE", reason=info["lead_guard_error"])
        return row, None
    N = dev.n_slices
    row["n_slices"], row["matrix_order"] = N, N * prob.Ny
    t0 = time.time()
    rg = scattering_matrix(E, dev)
    row["time_M1_s"] = time.time() - t0
    t0 = time.time()
    gb = sr.sparse_scattering_matrix(E, pot.numpy, geo["x_start"], dx, N, prob.Ly, prob.Ny, prob.Nq)
    row["time_M2_s"] = time.time() - t0
    t0 = time.time()
    wf = sw.wavefunction_scattering_matrix(E, pot.numpy, geo["x_start"], dx, N, prob.Ly, prob.Ny, prob.Nq)
    row["time_M3_s"] = time.time() - t0
    n = rg["n_open"]
    S = {"M1": rg["S"], "M2": gb["S"], "M3": wf["S"]}
    Ep = dev.basis.E_perp
    analytic = int(np.sum(0.5 * (np.arange(1, 5000) * np.pi / prob.Ly) ** 2 < E))
    row.update(n_open_M1=n, n_open_M2=gb["n_open"], n_open_M3=wf["n_open"], n_open_analytic=analytic, S_dim=int(rg["S"].shape[0]))
    row["N1_open_counts_equal"] = bool(n == gb["n_open"] == wf["n_open"] == analytic)
    for m, Sm in S.items():
        mm = b3x.method_metrics(Sm, n)
        row.update({f"{k}_{m}": v for k, v in mm.items()})
        row[f"T1_{m}"] = b3.T1(Sm, n)
        row[f"R1_{m}"] = float(np.sum(np.abs(b3.sm.split_S(Sm, n)["r_L"][:, 0]) ** 2))
    for (a, b, tag) in (("M1", "M2", "S1"), ("M1", "M3", "S2"), ("M2", "M3", "S3")):
        c = b3x.compare_S(S[a], S[b], n)
        row[f"{tag}_dS"] = c["dS"]
        row[f"{tag}_argmax"] = f"{c['argmax_row']},{c['argmax_col']}"
        if tag in ("S1", "S2"):
            row[f"S4_{tag}_dP_all"], row[f"S4_{tag}_dP_offdiag_trans"] = c["dP_all"], c["dP_offdiag_trans"]
            row[f"S4_{tag}_dT1"], row[f"S4_{tag}_dR1"] = c["dT1"], c["dR1"]
            row[f"S4_{tag}_max"] = max(c["dP_all"], c["dP_offdiag_trans"], c["dT1"], c["dR1"])
        row[f"rho_{tag}"] = rho_diag(c["dS"], E, Ep, dx, N, n)
    row["R1_backward_error_M2"] = gb["relative_residual"]
    row["R1_backward_error_M3"] = wf["residual_after_refinement"]
    row["M3_residual_before_refinement"] = wf["residual_before_refinement"]
    row["M1_min_rcond"] = rg["min_rcond"]
    cc = b3x.closed_channel_checks(E, rg, wf, Ep, dx)
    row.update({f"N1_{k}": v for k, v in cc.items()})
    row["N1_gamma_closed_M2"] = float(np.max(np.abs(gb["modes"]["gamma"][~gb["modes"]["open"]]))) if (~gb["modes"]["open"]).any() else 0.0
    row["N1_passed"] = bool(row["N1_open_counts_equal"] and cc["gamma_closed_max_rgf"] == 0.0 and cc["gamma_closed_max_wf"] == 0.0
                            and row["N1_gamma_closed_M2"] == 0.0 and cc["lam_closed_max_imag_rgf"] == 0.0 and cc["lam_closed_max_imag_wf"] == 0.0
                            and cc["lam_closed_max_abs_rgf"] < 1.0 and cc["lam_closed_max_abs_wf"] < 1.0 and cc["n_closed"] >= 1 and cc["open_masks_equal"])
    Ep2, Ep3 = sr.transverse_data(prob.Ly, prob.Ny, prob.Nq)[3], wf["E_perp"]
    row["P1_Eperp_bitwise_across_methods"] = bool(np.array_equal(Ep, Ep2) and np.array_equal(Ep, Ep3))
    row["P1_M2_samples_bitwise"] = bool(np.array_equal(gb["V_samples"], ref))
    row["P1_M3_samples_bitwise"] = bool(np.array_equal(wf["V_samples"], ref))
    row["V_block_rel_diff_M3_vs_M1"] = float(np.max(np.abs(wf["V_blocks"] - dev.V)) / np.max(np.abs(ref)))
    row["V_block_rel_diff_M3_vs_M2"] = float(np.max(np.abs(wf["V_blocks"] - gb["V_blocks"])) / np.max(np.abs(ref)))
    row["status"] = "OK"
    return row, (S if keep_S else None)


# ======================================================================
# Tier 3 (M4) for one point and resolution, with the decomposition for M2 and M3
# ======================================================================
def run_xp_case(prob, devs, pt, inv, buffer, S_double, decompose=False, log=print):
    cid, E = pt["config"], pt["E"]
    dev, geo, ref, info, pot = devs.get(cid, inv, buffer)
    dx = 1.0 / inv
    Ep = dev.basis.E_perp
    t0 = time.time()
    ex = xp.exact_lattice_S(E, Ep, ref, dx, prob.Ny, prob.Nq)
    S_ex = ex["S"]
    row = dict(point=pt["id"], config=cid, E=E, inv_dx=inv, dx=dx, buffer_requested=buffer, n_slices=ref.shape[0], n_open=ex["n_open"],
               potential_hash=info["potential_hash"], time_M4_s=time.time() - t0,
               X4_state_passed=ex["refinement_gate"]["passed"], X4_state_worst_contraction=ex["refinement_gate"]["worst_contraction"],
               X4_state_iterations=ex["refinement_gate"]["n_iterations"], X4_state_eta_final=ex["eta_final"],
               X4_green_passed=ex["green_refinement_gate"]["passed"], X4_green_worst_contraction=ex["green_refinement_gate"]["worst_contraction"],
               X4_green_eta_final=ex["green_refinement_gate"]["eta_final"], X5_fl_vs_state_max_abs=ex["fl_vs_state_max_abs"],
               X5_passed=bool(ex["fl_vs_state_max_abs"] <= xp.FL_WF_TOL), X7_error_indicator=ex["error_indicator"],
               X7_passed=bool(ex["error_indicator"] <= xp.ERROR_INDICATOR_MAX), kappa1=ex["kappa1"],
               T1_M4=float(np.sum(np.abs(ex["t_RL"][:, 0]) ** 2)))
    row["gates_passed"] = bool(row["X4_state_passed"] and row["X4_green_passed"] and row["X5_passed"] and row["X7_passed"])
    for m, Sm in S_double.items():
        d = np.abs(Sm.astype(CLD) - S_ex)
        i = np.unravel_index(int(np.argmax(d)), d.shape)
        row[f"eps_{m}"] = float(d.max())
        row[f"eps_{m}_argmax"] = f"{i[0]},{i[1]}"
        row[f"dT1_{m}_vs_M4"] = abs(b3.T1(Sm, ex["n_open"]) - row["T1_M4"])
    out = dict(row=row, S_exact=S_ex, decomp=[])
    if decompose and row["gates_passed"]:
        base = dict(point=pt["id"], config=cid, E=E, inv_dx=inv, dx=dx, n_slices=ref.shape[0], n_open=ex["n_open"])
        t1 = time.time()
        wf = sw.wavefunction_scattering_matrix(E, pot.numpy, devs.get(cid, inv, buffer)[1]["x_start"], dx, ref.shape[0], prob.Ly, prob.Ny, prob.Nq)
        gb = sr.sparse_scattering_matrix(E, pot.numpy, devs.get(cid, inv, buffer)[1]["x_start"], dx, ref.shape[0], prob.Ly, prob.Ny, prob.Nq)
        rep = dict(M2_recomputed_bitwise=bool(np.array_equal(gb["S"], S_double["M2"])), M3_recomputed_bitwise=bool(np.array_equal(wf["S"], S_double["M3"])))
        for d in (decompose_m2(E, gb, sr.transverse_data(prob.Ly, prob.Ny, prob.Nq)[3], dx, ex["lead"], S_ex, ex["error_indicator"]),
                  decompose_m3(E, wf, dx, ex["lead"], S_ex, ex["error_indicator"])):
            r = decomposition_row(dict(base, **rep), d)
            r["time_decomposition_s"] = time.time() - t1
            out["decomp"].append(r)
    elif decompose:
        for m in ("M2", "M3"):
            out["decomp"].append(dict(point=pt["id"], inv_dx=inv, method=m, status="NOT EVALUABLE", reason="extended-precision gates X4/X5/X7 failed at this point"))
    return out


# ======================================================================
# gates X0-X3, X6 (extended precision vs mpmath / stability), G0, potential identity P2
# ======================================================================
def gate_X0():
    p = xp.platform_precision()
    return p, bool(p["evaluable"])


def eperp_of(prob):
    return 0.5 * (np.arange(1, prob.Ny + 1) * np.pi / prob.Ly) ** 2


def gate_X1(prob, points, log=print):
    rows = []
    Ep = eperp_of(prob)
    for pt in points:
        for inv in prob.invs:
            r = xp.lead_data_vs_mpmath(pt["E"], Ep, 1.0 / inv)
            rows.append(dict(gate="X1", case=pt["id"], inv_dx=inv, value=max(r["rel_lambda"], r["rel_gamma"]), rel_lambda=r["rel_lambda"],
                             rel_gamma=r["rel_gamma"], limit=XPG["lead_data_vs_mpmath_rel"], passed=r["passed"]))
    return rows


def _clean_or_barrier_cases(prob, pts, invs, kind, log=print):
    """X2 (clean guide) and X3 (y-uniform smooth barrier) on the geometry of the predeclared devices at the given points' energies."""
    import mpmath as mp
    rows = []
    Ep = eperp_of(prob)
    cid0 = next(iter(pts))["config"]
    for pt in pts:
        for inv in invs:
            geo = prob.geometry(pt["config"], inv, prob.buffer_req)
            dx, N = 1.0 / inv, geo["n_slices"]
            xs = geo["x_start"] + dx * np.arange(N)
            if kind == "X2":
                Vx = np.zeros(N)
            else:
                a, b = prob.region(pt["config"])
                c, w = 0.5 * (a + b), 0.5 * (b - a) * 0.1
                Vx = 150.0 * 0.5 * (np.tanh((xs - (c - w)) / 0.01) - np.tanh((xs - (c + w)) / 0.01))
            V = np.tile(Vx[:, None], (1, prob.Nq))
            t0 = time.time()
            res = xp.exact_lattice_S(pt["E"], Ep, V, dx, prob.Ny, prob.Nq)
            no = res["n_open"]
            worst_amp, worst_off = 0.0, 0.0
            with mp.workdps(xp.MPMATH_DPS):
                if kind == "X2":
                    ref = xp.lead_data_mpmath(pt["E"], Ep, dx)
                    for n in range(no):
                        ex = ref["lam"][n] ** (N - 1)
                        worst_amp = max(worst_amp, float(abs(xp._to_mp(res["t_RL"][n, n]) - ex)), float(abs(xp._to_mp(res["t_LR"][n, n]) - ex)))
                    worst_off = float(max(np.max(np.abs(res["t_RL"] - np.diag(np.diag(res["t_RL"])))), np.max(np.abs(res["r_L"])), np.max(np.abs(res["r_R"]))))
                else:
                    for n in range(no):
                        r_ex, t_ex = xp.mp_single_channel_amplitudes(pt["E"], Ep[n], Vx, dx)
                        worst_amp = max(worst_amp, float(abs(xp._to_mp(res["t_RL"][n, n]) - t_ex)), float(abs(xp._to_mp(res["r_L"][n, n]) - r_ex)))
                    worst_off = float(np.max(np.abs(res["t_RL"] - np.diag(np.diag(res["t_RL"])))))
            g = res["refinement_gate"]
            gates = bool(g["passed"] and res["green_refinement_gate"]["passed"] and res["fl_vs_state_max_abs"] <= xp.FL_WF_TOL
                         and res["error_indicator"] <= xp.ERROR_INDICATOR_MAX)
            rows.append(dict(gate=kind, case=pt["id"], inv_dx=inv, value=max(worst_amp, worst_off), amplitude_error_vs_mpmath=worst_amp,
                             channel_mixing_or_reflection=worst_off, limit=XPG["scalar_problem_vs_mpmath_abs"],
                             passed=bool(max(worst_amp, worst_off) <= XPG["scalar_problem_vs_mpmath_abs"] and gates), X4_X5_X7_passed=gates,
                             eta_final=res["eta_final"], e_ind=res["error_indicator"], time_s=time.time() - t0))
            log(f"    {kind} {pt['id']} 1/{inv}: {rows[-1]['value']:.2e} gates={gates} ({rows[-1]['time_s']:.0f}s)")
    return rows


def gate_X6(prob, devs, pts, invs, log=print):
    rows = []
    Ep = eperp_of(prob)
    for pt in pts:
        for inv in invs:
            dev, geo, ref, info, pot = devs.get(pt["config"], inv, prob.buffer_req)
            Vb = xp.potential_blocks_extended(ref, prob.Ny, prob.Nq)
            r = xp.ordering_stability(pt["E"], Ep, Vb, 1.0 / inv)
            rows.append(dict(gate="X6", case=pt["id"], inv_dx=inv, value=r["max_abs_difference"], limit=xp.ORDERING_TOL, passed=r["passed"],
                             orderings="/".join(r["orderings"])))
            log(f"    X6 {pt['id']} 1/{inv}: {r['max_abs_difference']:.2e}")
    return rows


def potential_identity_rows(prob, devs, pts_by_inv):
    rows = []
    for (pid, cid, inv, buf) in pts_by_inv:
        dev, geo, ref, info, pot = devs.get(cid, inv, buf)
        rows.append(dict(config=cid, inv_dx=inv, buffer_requested=buf, n_slices=ref.shape[0], x_start=geo["x_start"], x_end=geo["x_end"], shape=str(info["shape"]),
                         potential_hash=info["potential_hash"], x_hash=info["x_hash"], V_min=info["V_min"], V_max=info["V_max"], V_lead_max=info["V_lead_max"],
                         lead_guard_limit=LIMITS["lead_tail_guard"], P4_lead_guard_ok=bool(info["V_lead_max"] <= LIMITS["lead_tail_guard"]),
                         P1_rgf_device_vs_direct_bitwise=info["bitwise_equal"], P1_max_abs_diff=info["max_abs_diff"]))
    return rows


def cross_resolution_identity(prob, devs, cids, invs):
    """P2: V at coincident physical slices is bitwise equal across resolutions (same buffer length, nested lattices)."""
    rows = []
    for cid in cids:
        base_inv = invs[0]
        _, geo0, ref0, _, _ = devs.get(cid, base_inv, prob.buffer_req)
        x0 = geo0["x_start"] + (1.0 / base_inv) * np.arange(ref0.shape[0])
        for inv in invs[1:]:
            _, geo, ref, _, _ = devs.get(cid, inv, prob.buffer_req)
            step = inv // base_inv
            sub = ref[::step][:ref0.shape[0]]
            xs = geo["x_start"] + (1.0 / inv) * np.arange(ref.shape[0])
            same_x = bool(np.array_equal(xs[::step][:ref0.shape[0]], x0))
            rows.append(dict(config=cid, inv_dx_fine=inv, inv_dx_coarse=base_inv, n_coincident=int(ref0.shape[0]), coincident_x_bitwise=same_x,
                             P2_potential_bitwise=bool(same_x and sub.shape == ref0.shape and np.array_equal(sub, ref0)),
                             max_abs_diff=float(np.max(np.abs(sub - ref0))) if sub.shape == ref0.shape else float("nan")))
    return rows


# ======================================================================
# caching (results are assembled from per-case caches so a restart loses at most one case)
# ======================================================================
class Cache:
    def __init__(self, path, stamp):
        self.path = path
        os.makedirs(path, exist_ok=True)
        sp = os.path.join(path, "cache_stamp.json")
        if os.path.exists(sp):
            with open(sp) as f:
                old = json.load(f)
            if old != stamp:
                raise RuntimeError(f"cache {path} was written for a different protocol/problem ({old} vs {stamp}); use a fresh --cache directory")
        else:
            with open(sp, "w") as f:
                json.dump(stamp, f)

    def get(self, key, fn):
        p = os.path.join(self.path, key.replace("/", "_") + ".pkl")
        if os.path.exists(p):
            with open(p, "rb") as f:
                return pickle.load(f), True
        v = fn()
        tmp = p + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(v, f)
        os.replace(tmp, p)
        return v, False


# ======================================================================
# campaign
# ======================================================================
def check_output_dir(out):
    o = os.path.normpath(os.path.abspath(out))
    for d in FORBIDDEN_OUTPUT_DIRS:
        if o == d or o.startswith(d + os.sep):
            raise SystemExit(f"refusing to write into the historical result directory {d}")
    return o


def load_b3x_rows():
    with open(B3X_CSV) as f:
        return list(csv.DictReader(f))


def stage_G0(rows512, log=print):
    """G0: the recomputed 1/512 values reproduce the B3X artifacts (RGF T1 bitwise; rgf_wf / rgf_gb / wf_gb dS within 1e-13)."""
    b3 = load_b3x_rows()
    out = []
    for pid, (row, S) in rows512.items():
        cand = [r for r in b3 if r["config"] == row["config"] and float(r["E"]) == row["E"]]
        if len(cand) != 1:
            out.append(dict(gate="G0", case=pid, passed=False, value=float("nan"), reason=f"{len(cand)} matching B3X rows"))
            continue
        r = cand[0]
        t1_bit = float(r["T1_rgf"]) == row["T1_M1"]
        d = [abs(float(r["rgf_wf_dS"]) - row["S2_dS"]), abs(float(r["rgf_gb_dS"]) - row["S1_dS"]), abs(float(r["wf_gb_dS"]) - row["S3_dS"])]
        out.append(dict(gate="G0", case=pid, value=max(d), T1_rgf_bitwise=t1_bit, d_rgf_wf=d[0], d_rgf_gb=d[1], d_wf_gb=d[2], limit=1e-13,
                        passed=bool(t1_bit and max(d) <= 1e-13)))
    return out


def run_campaign(prob, out, cache_dir, stages=None, log=print, run_g0=True):
    out = check_output_dir(out)
    os.makedirs(out, exist_ok=True)
    stamp = dict(problem=prob.name, protocol_sha256=sha_file(PROTOCOL_PATH), invs=list(prob.invs))
    cache = Cache(cache_dir, stamp)
    devs = DeviceCache(prob)
    t_start = time.time()
    P = prob.points
    sets = prob.sets_map
    R = dict(comparison=[], xp=[], decomp=[], gates=[], identity=[], p2=[], buffer=[], notes=[], stopped_at=None, xp_available=True)
    rows = {}                                   # (pid, inv) -> (row, S)

    def methods(pid, inv, buf=None):
        buf = prob.buffer_req if buf is None else buf
        key = f"M123__{pid}__{inv}__{buf}"
        (row, S), hit = cache.get(key, lambda: run_methods(prob, devs, P[pid], inv, buf))
        return row, S

    def stop(stage, why):
        R["stopped_at"] = dict(stage=stage, reason=why)
        log(f"STOP at {stage}: {why}")

    # ---- G0 + Tier 1 at the base resolution (also the 1/512 validation inputs) ----
    log("stage G0 / Tier 1 @ base resolution")
    inv0 = prob.inv_base
    base_ids = list(dict.fromkeys(sets["tier1"] + sets["tier3_validation"]))
    for pid in base_ids:
        rows[(pid, inv0)] = methods(pid, inv0)
    R["comparison"] += [rows[(pid, inv0)][0] for pid in base_ids]
    if run_g0:
        g0 = stage_G0({pid: rows[(pid, inv0)] for pid in base_ids})
        R["gates"] += g0
        if not all(g["passed"] for g in g0):
            stop("G0", "recomputed 1/512 values do not reproduce the B3X artifacts")
            return finish(prob, out, R, rows, t_start, cache_dir, log)
    # ---- potential identity ----
    log("stage potential identity (P1-P4)")
    cids = sorted({P[p]["config"] for p in base_ids})
    ident_keys = [(p, P[p]["config"], inv0, prob.buffer_req) for p in base_ids]
    # ---- Tier 1 at the finer levels ----
    mid, fine = prob.invs[1], prob.invs[2]
    log(f"stage Tier 1 @ 1/{mid}")
    for pid in sets["tier1"]:
        rows[(pid, mid)] = methods(pid, mid)
        R["comparison"].append(rows[(pid, mid)][0])
    log(f"stage Tier 2 @ 1/{fine}")
    for pid in sets["tier2"]:
        rows[(pid, fine)] = methods(pid, fine)
        R["comparison"].append(rows[(pid, fine)][0])
    # 1/512 rows of the five points also enter the convergence study (already computed)
    ident_keys += [(p, P[p]["config"], mid, prob.buffer_req) for p in sets["tier1"]] + [(p, P[p]["config"], fine, prob.buffer_req) for p in sets["tier2"]]
    R["identity"] = potential_identity_rows(prob, devs, list(dict.fromkeys(ident_keys)))
    R["p2"] = cross_resolution_identity(prob, devs, sorted({P[p]["config"] for p in sets["tier2"]}), prob.invs)
    # ---- XP gates ----
    log("stage XP gates (X0-X3, X6)")
    p, ok0 = gate_X0()
    R["platform"] = p
    R["gates"].append(dict(gate="X0", case="platform", value=p.get("longdouble_mantissa_bits"), limit=XPG["platform_min_mantissa_bits"], passed=ok0))
    if not ok0:
        R["xp_available"] = False
        stop("X0", "NOT EVALUABLE: extended precision unavailable; Tier 3 not run (no fallback to float64)")
    else:
        tier3_pts = [P[p_] for p_ in dict.fromkeys(sets["tier3_validation"] + sets["tier3_fine"])]
        fivepts = [P[p_] for p_ in sets["tier3_fine"]]
        g, _ = cache.get("gate_X1", lambda: gate_X1(prob, tier3_pts, log))
        R["gates"] += g
        for kind in ("X2", "X3"):
            g, _ = cache.get(f"gate_{kind}", lambda kind=kind: _clean_or_barrier_cases(prob, fivepts, prob.invs, kind, log))
            R["gates"] += g
        g, _ = cache.get("gate_X6", lambda: gate_X6(prob, devs, fivepts, prob.invs[:2], log))
        R["gates"] += g
        failed = [x for x in R["gates"] if x["gate"] in ("X0", "X1", "X2", "X3", "X6") and not x["passed"]]
        if failed:
            R["xp_available"] = False
            stop("XP gates", f"{len(failed)} extended-precision gate rows failed ({sorted({x['gate'] for x in failed})}); XP-dependent stages skipped")
    # ---- X8 / V1 and Tier 3 ----
    xp_rows = {}
    if R["xp_available"]:
        log("stage X8 / V1: M4 vs M1-M3 at the base resolution")
        for pid in sets["tier3_validation"]:
            deco = pid in sets["tier3_fine"]
            res, _ = cache.get(f"M4__{pid}__{inv0}", lambda pid=pid, deco=deco: run_xp_case(prob, devs, P[pid], inv0, prob.buffer_req, rows[(pid, inv0)][1], decompose=deco, log=log))
            xp_rows[(pid, inv0)] = res
            log(f"  {pid}: eps M1/M2/M3 = {res['row']['eps_M1']:.2e} / {res['row']['eps_M2']:.2e} / {res['row']['eps_M3']:.2e}  gates={res['row']['gates_passed']}")
        v1 = [xp_rows[(pid, inv0)]["row"] for pid in sets["tier3_validation"]]
        x8 = all(r["gates_passed"] and max(r["eps_M1"], r["eps_M2"], r["eps_M3"]) <= XPG["validation_vs_B3X_refs_at_dx_1_512_abs"] for r in v1)
        R["gates"].append(dict(gate="X8", case="all validation points @ base resolution", value=max(max(r["eps_M1"], r["eps_M2"], r["eps_M3"]) for r in v1),
                               limit=XPG["validation_vs_B3X_refs_at_dx_1_512_abs"], passed=bool(x8)))
        if not x8:
            R["xp_available"] = False
            stop("X8/V1", "Tier 3 not validated at the base resolution; XP is not used to interpret finer grids")
    if R["xp_available"]:
        for inv in (mid, fine):
            log(f"stage Tier 3 @ 1/{inv}")
            for pid in sets["tier3_fine"]:
                res, _ = cache.get(f"M4__{pid}__{inv}", lambda pid=pid, inv=inv: run_xp_case(prob, devs, P[pid], inv, prob.buffer_req, rows[(pid, inv)][1] if (pid, inv) in rows else
                                                                                           methods(pid, inv)[1], decompose=True, log=log))
                xp_rows[(pid, inv)] = res
                r = res["row"]
                log(f"  {pid} 1/{inv}: eps M1/M2/M3 = {r['eps_M1']:.2e} / {r['eps_M2']:.2e} / {r['eps_M3']:.2e}  gates={r['gates_passed']}  ({r['time_M4_s']:.0f}s)")
    for k, res in xp_rows.items():
        R["xp"].append(res["row"])
        R["decomp"] += res["decomp"]
    # ---- buffer control ----
    log("stage buffer control K1")
    cbuf = prob.buffer_ctrl_req
    ibuf = []
    for pid in sets["buffer"]:
        rows[(pid, mid, cbuf)] = methods(pid, mid, cbuf)
        ibuf.append((pid, P[pid]["config"], mid, cbuf))
    R["identity"] += potential_identity_rows(prob, devs, ibuf)
    for pid in sets["buffer"]:
        r103, S103 = rows[(pid, mid)]
        rc, Sc = rows[(pid, mid, cbuf)]
        if r103.get("status") != "OK" or rc.get("status") != "OK":
            R["buffer"].append(dict(point=pid, status="NOT EVALUABLE", reason="lead-tail guard"))
            continue
        row = dict(point=pid, config=P[pid]["config"], E=P[pid]["E"], inv_dx=mid, buffer_default=prob.buffer_req, buffer_control=cbuf,
                   n_slices_default=r103["n_slices"], n_slices_control=rc["n_slices"], buffer_length_default=r103["buffer_length"], buffer_length_control=rc["buffer_length"],
                   limit_dT1=LIMITS["dT1_buffer"], limit_S=LIMITS["S_abs"])
        for m in ("M1", "M3"):
            row[f"T1_{m}_default"], row[f"T1_{m}_control"] = r103[f"T1_{m}"], rc[f"T1_{m}"]
            row[f"abs_dT1_{m}"] = abs(rc[f"T1_{m}"] - r103[f"T1_{m}"])
            row[f"dT1_{m}_passed"] = passes(row[f"abs_dT1_{m}"], LIMITS["dT1_buffer"])
        row["S_M1_vs_M3_default"], row["S_M1_vs_M3_control"] = r103["S2_dS"], rc["S2_dS"]
        row["S_within_geometry_passed"] = bool(passes(r103["S2_dS"], LIMITS["S_abs"]) and passes(rc["S2_dS"], LIMITS["S_abs"]))
        row["abs_d_probabilities_M1"] = float(np.max(np.abs(np.abs(S103["M1"]) ** 2 - np.abs(Sc["M1"]) ** 2)))
        row["K1_passed_dT1"] = bool(row["dT1_M1_passed"] and row["dT1_M3_passed"])
        R["buffer"].append(row)
    return finish(prob, out, R, rows, t_start, cache_dir, log, xp_rows)


# ======================================================================
# analysis, summary and output
# ======================================================================
def _max_arg(rows, key):
    v = [(r[key], r["point"], r["inv_dx"]) for r in rows if key in r and r.get(key) is not None and np.isfinite(r[key])]
    return max(v) if v else (float("nan"), None, None)


def row_status(r):
    """Per-row criterion statuses (predeclared vocabulary)."""
    if r.get("status") != "OK":
        return {k: "NOT EVALUABLE" for k in ("S1", "S2", "S3", "S4", "S5", "S6", "R1", "N1")}
    st = {}
    for tag in ("S1", "S2", "S3"):
        st[tag] = "PASS" if passes(r[f"{tag}_dS"], LIMITS["S_abs"]) else "PRECISION-LIMITED (inconclusive)"
    s4 = max(r["S4_S1_max"], r["S4_S2_max"])
    if passes(s4, LIMITS["TR_offdiag_abs"]):
        st["S4"] = "PASS"
    else:
        st["S4"] = "PRECISION-LIMITED (inconclusive)" if (st["S1"] != "PASS" or st["S2"] != "PASS") else "FAIL"
    st["S5"] = "PASS" if all(passes(r[f"U_{m}"], LIMITS["unitarity"]) for m in ("M1", "M2", "M3")) else "FAIL"
    st["S6"] = "PASS" if all(passes(r[f"Q_{m}"], LIMITS["reciprocity"]) for m in ("M1", "M2", "M3")) else "FAIL"
    st["R1"] = "PASS" if passes(max(r["R1_backward_error_M2"], r["R1_backward_error_M3"]), LIMITS["sparse_backward_error"]) else "FAIL"
    st["N1"] = "PASS" if r["N1_passed"] else "FAIL"
    return st


def convergence_rows(prob, rows, xp_rows):
    """O1: T1 at 1/512, 1/1024, 1/2048 per method for the Tier 2 points. The Richardson value is an ESTIMATE, applied only if the observed order is in the window."""
    out = []
    invs = prob.invs
    lo, hi = LIMITS["order_window"]
    for pid in prob.sets_map["tier2"]:
        for m in ("M1", "M2", "M3", "M4"):
            try:
                if m == "M4":
                    T = [xp_rows[(pid, i)]["row"]["T1_M4"] for i in invs]
                else:
                    T = [rows[(pid, i)][0][f"T1_{m}"] for i in invs]
            except KeyError:
                continue
            d1, d2 = T[1] - T[0], T[2] - T[1]
            p_obs = math.log2(abs(d1) / abs(d2)) if d1 != 0 and d2 != 0 else float("nan")
            ok = bool(np.isfinite(p_obs) and lo <= p_obs <= hi)
            est = T[2] + d2 / (2 ** p_obs - 1) if ok else float("nan")
            out.append(dict(point=pid, method=m, **{f"T1_{i}": t for i, t in zip(invs, T)}, delta_1=d1, delta_2=d2, observed_order=p_obs,
                            n_orders_available=1, order_in_window=ok, window=f"[{lo},{hi}]", T1_richardson_ESTIMATE=est,
                            note="only ONE observed order exists with three grids; the protocol's 'both orders' rule is applied to it (reported ambiguity); estimate, not a continuum claim"))
    return out


def solver_spread(prob, rows, xp_rows):
    out = []
    for pid in prob.sets_map["tier2"]:
        for inv in prob.invs:
            if (pid, inv) not in rows:
                continue
            r = rows[(pid, inv)][0]
            T = [r["T1_M1"], r["T1_M2"], r["T1_M3"]]
            if (pid, inv) in xp_rows:
                T.append(xp_rows[(pid, inv)]["row"]["T1_M4"])
            out.append(dict(point=pid, inv_dx=inv, n_methods=len(T), T1_spread=max(T) - min(T), limit=LIMITS["T1_method_spread_for_solver_independence"],
                            passed=passes(max(T) - min(T), LIMITS["T1_method_spread_for_solver_independence"])))
    return out


def adjudication(prob, xp_rows):
    """Tier 3 adjudication (protocol step 1-6) for the Tier 2 / Tier 3 fine points."""
    out = []
    for pid in prob.sets_map["tier3_fine"]:
        for m in ("M1", "M2", "M3"):
            canon = {inv: c for inv, c in zip(prob.invs, (512, 1024, 2048))}       # classification rules are stated for 1/512, 1/1024, 1/2048 (identity for B3Y)
            eps = {}
            for inv in prob.invs:
                if (pid, inv) in xp_rows and xp_rows[(pid, inv)]["row"]["gates_passed"]:
                    eps[canon[inv]] = xp_rows[(pid, inv)]["row"][f"eps_{m}"]
            if not eps:
                continue
            back = {c: i for i, c in canon.items()}
            exceed = {i: v for i, v in eps.items() if v > LIMITS["S_abs"]}
            pc = xp.precision_consistency(eps) if len(eps) >= 2 else dict(eps_512_within_target=False, growth_per_halving=[], passed=False)
            entry = dict(point=pid, method=m, eps_by_inv_dx=eps, exceeds_S_target_at=sorted(exceed), precision_consistency=pc)
            if not exceed:
                entry["classification"] = "CONSISTENT WITH THE EXACT LATTICE SOLUTION AT 1e-10 (all available resolutions)"
            else:
                label = xp.classify_exceedance(eps, None, float("nan"), claim_representation_dominance=False)
                entry["classification"] = label
                entry["mechanism_claim"] = "none"
                if m in ("M2", "M3"):
                    fin = max(exceed)
                    dr = [d for d in xp_rows[(pid, back[fin])]["decomp"] if d.get("method") == m and "bulk_rounding" in d]
                    if dr:
                        d = dr[0]
                        entry["decomposition_at_inv_dx"] = fin
                        entry["decomposition_resolved"] = d["resolved"]
                        entry["decomposition_dominance_supported"] = d["dominance_supported"]
                        entry["decomposition_reason"] = d["dominance_reason"]
                        if label == "PRECISION-CONSISTENT":
                            entry["classification"] = ("PRECISION-LIMITED (measured; representation-dominated)" if d["dominance_supported"]
                                                       else "PRECISION-CONSISTENT (decomposition does not support a representation-dominance claim)")
                            entry["mechanism_claim"] = "representation dominance supported by the decomposition" if d["dominance_supported"] else "none"
                    else:
                        entry["decomposition_reason"] = "no decomposition available at this point/resolution"
                else:
                    entry["decomposition_reason"] = "M1 has no explicit matrix: only its total error is available; no mechanism claim"
            out.append(entry)
    return out


def finish(prob, out, R, rows, t_start, cache_dir, log, xp_rows=None):
    xp_rows = xp_rows or {}
    comp = R["comparison"]
    for r in comp:
        r.update({f"status_{k}": v for k, v in row_status(r).items()})
    conv = convergence_rows(prob, rows, xp_rows) if all((p, i) in rows for p in prob.sets_map["tier2"] for i in prob.invs) else []
    spread = solver_spread(prob, rows, xp_rows)
    adj = adjudication(prob, xp_rows)
    write_csv(os.path.join(out, "per_energy_comparison.csv"), comp)
    write_csv(os.path.join(out, "tier3_exact_lattice.csv"), R["xp"])
    write_csv(os.path.join(out, "error_decomposition.csv"), R["decomp"])
    write_csv(os.path.join(out, "refinement_convergence.csv"), conv)
    write_csv(os.path.join(out, "buffer_control.csv"), R["buffer"])
    write_csv(os.path.join(out, "potential_identity.csv"), R["identity"] + [dict(config=x["config"], inv_dx=f"{x['inv_dx_fine']}_vs_{x['inv_dx_coarse']}", P2_potential_bitwise=x["P2_potential_bitwise"],
                                                                                coincident_x_bitwise=x["coincident_x_bitwise"], max_abs_diff=x["max_abs_diff"], n_slices=x["n_coincident"]) for x in R["p2"]])
    write_csv(os.path.join(out, "xp_gates.csv"), R["gates"])
    write_csv(os.path.join(out, "solver_spread.csv"), spread)
    summ = summarize(prob, R, comp, conv, spread, adj, xp_rows)
    summ["stopped_at"] = R["stopped_at"]
    summ["wall_time_s"] = time.time() - t_start
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summ, f, indent=2, default=_json_default)
    with open(os.path.join(out, "platform_precision.json"), "w") as f:
        json.dump(platform_precision_report(R.get("platform")), f, indent=2, default=_json_default)
    write_metadata(prob, out, R, comp, summ, cache_dir)
    log(f"wrote results to {out} ({summ['wall_time_s']:.0f}s)")
    return summ


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def summarize(prob, R, comp, conv, spread, adj, xp_rows):
    ok = [r for r in comp if r.get("status") == "OK"]
    s = dict(problem=prob.name, n_comparison_rows=len(comp), n_not_evaluable=len(comp) - len(ok), protocol_status=PROTOCOL["status"], criteria={})
    crit = {}
    for tag, key, lim in (("S1_M1_vs_M2", "S1_dS", LIMITS["S_abs"]), ("S2_M1_vs_M3", "S2_dS", LIMITS["S_abs"]), ("S3_M2_vs_M3", "S3_dS", LIMITS["S_abs"]),
                          ("S4_TR_offdiag", None, LIMITS["TR_offdiag_abs"])):
        if key is None:
            vals = [(max(r["S4_S1_max"], r["S4_S2_max"]), r) for r in ok]
        else:
            vals = [(r[key], r) for r in ok]
        v, r = max(vals, key=lambda x: x[0]) if vals else (float("nan"), None)
        nexc = sum(1 for x, _ in vals if not passes(x, lim))
        crit[tag] = dict(limit=lim, max=v, argmax=(r["point"], r["inv_dx"]) if r else None, n=len(vals), n_exceed=nexc,
                         status="PASS" if nexc == 0 else "PRECISION-LIMITED (inconclusive)" if tag != "S4_TR_offdiag" or any(
                             row_status(x)["S4"] != "FAIL" for _, x in vals if not passes(_, lim)) else "FAIL")
    for tag, st_key in (("S5_unitarity", "S5"), ("S6_reciprocity", "S6"), ("R1_backward_error", "R1"), ("N1_channels", "N1")):
        sts = [row_status(r)[st_key] for r in ok]
        crit[tag] = dict(n=len(sts), n_fail=sum(x == "FAIL" for x in sts), status="PASS" if all(x == "PASS" for x in sts) else "FAIL")
    crit["S5_unitarity"]["max"] = max((max(r[f"U_{m}"] for m in ("M1", "M2", "M3")) for r in ok), default=float("nan"))
    crit["S6_reciprocity"]["max"] = max((max(r[f"Q_{m}"] for m in ("M1", "M2", "M3")) for r in ok), default=float("nan"))
    crit["R1_backward_error"]["max"] = max((max(r["R1_backward_error_M2"], r["R1_backward_error_M3"]) for r in ok), default=float("nan"))
    p1 = [all([r["P1_rgf_device_bitwise"], r["P1_M2_samples_bitwise"], r["P1_M3_samples_bitwise"], r["P1_Eperp_bitwise_across_methods"]]) for r in ok]
    crit["P1_potential_identity"] = dict(n=len(p1), status="PASS" if all(p1) and all(x["P1_rgf_device_vs_direct_bitwise"] for x in R["identity"]) else "FAIL")
    crit["P2_cross_resolution"] = dict(n=len(R["p2"]), status="PASS" if R["p2"] and all(x["P2_potential_bitwise"] for x in R["p2"]) else "FAIL")
    crit["P4_lead_guard"] = dict(n=len(comp), n_violating=sum(not r.get("P4_lead_guard_ok", False) for r in comp),
                                 status="PASS" if all(r.get("P4_lead_guard_ok", False) for r in comp) else "NOT EVALUABLE")
    g0 = [g for g in R["gates"] if g["gate"] == "G0"]
    crit["G0_reproduces_B3X"] = dict(n=len(g0), status="PASS" if g0 and all(g["passed"] for g in g0) else ("not run (synthetic)" if not g0 else "FAIL"))
    gates = {}
    for g in R["gates"]:
        gates.setdefault(g["gate"], []).append(g["passed"])
    crit["XP_gates"] = {k: dict(n=len(v), status="PASS" if all(v) else "FAIL") for k, v in gates.items()}
    crit["V1_tier3_validated"] = dict(status="PASS" if gates.get("X8") and all(gates["X8"]) and all(all(v) for k, v in gates.items() if k in ("X0", "X1", "X2", "X3", "X6")) else
                                      ("NOT EVALUABLE" if R["stopped_at"] and R["stopped_at"]["stage"] == "X0" else "FAIL"))
    if R["buffer"]:
        b = [x for x in R["buffer"] if "abs_dT1_M1" in x]
        crit["K1_buffer_control"] = dict(n=len(b), max_abs_dT1=max((max(x["abs_dT1_M1"], x["abs_dT1_M3"]) for x in b), default=float("nan")), limit=LIMITS["dT1_buffer"],
                                         status="PASS" if b and all(x["K1_passed_dT1"] and x["S_within_geometry_passed"] for x in b) else "FAIL")
    crit["H1_solver_independence_T1_spread"] = dict(n=len(spread), max=max((x["T1_spread"] for x in spread), default=float("nan")), limit=LIMITS["T1_method_spread_for_solver_independence"],
                                                    status="PASS" if spread and all(x["passed"] for x in spread) else "FAIL" if spread else "no data")
    crit["O1_convergence_orders"] = dict(status="CLASSIFICATION ONLY", orders={f"{c['point']}/{c['method']}": c["observed_order"] for c in conv},
                                         in_window={f"{c['point']}/{c['method']}": c["order_in_window"] for c in conv})
    s["criteria"] = crit
    s["tier3_adjudication"] = adj
    exc = [k for k, v in crit.items() if isinstance(v, dict) and v.get("status") in ("PRECISION-LIMITED (inconclusive)",)]
    hard_fail = [k for k, v in crit.items() if isinstance(v, dict) and v.get("status") == "FAIL"]
    nev = [k for k, v in crit.items() if isinstance(v, dict) and v.get("status") == "NOT EVALUABLE"]
    if hard_fail:
        overall = "FAIL"
    elif nev:
        overall = "NOT EVALUABLE (partial)"
    elif exc:
        if adj and all(a["classification"].startswith(("CONSISTENT", "PRECISION-CONSISTENT", "PRECISION-LIMITED")) for a in adj):
            overall = "PRECISION-LIMITED / PRECISION-CONSISTENT (resolved by Tier 3 where available); not PASS"
        else:
            overall = "INCONCLUSIVE / UNEXPLAINED (exceedances not fully adjudicated); not PASS"
    else:
        overall = "PASS"
    s["overall_status"] = overall
    s["overall_status_rule"] = PROTOCOL["tier3_adjudication"]["overall"]
    s["exceedances_1e-10"] = [dict(point=r["point"], inv_dx=r["inv_dx"], **{k: r[k] for k in ("S1_dS", "S2_dS", "S3_dS")}) for r in ok
                              if any(not passes(r[k], LIMITS["S_abs"]) for k in ("S1_dS", "S2_dS", "S3_dS"))]
    return s


def sha_file(path):
    with open(path, "rb") as f:
        return "sha256:" + hashlib.sha256(f.read()).hexdigest()


def git_out(*args):
    try:
        return subprocess.run(["git", *args], capture_output=True, text=True, cwd=HERE, timeout=20).stdout.strip() or None
    except Exception:                                   # noqa: BLE001
        return None


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


def platform_precision_report(prec):
    import mpmath
    f64, fl = np.finfo(np.float64), np.finfo(np.longdouble)
    cpu = None
    try:
        with open("/proc/cpuinfo") as f:
            cpu = next((l.split(":", 1)[1].strip() for l in f if l.startswith("model name")), None)
    except OSError:
        pass
    mem = None
    try:
        with open("/proc/meminfo") as f:
            mem = next((l.split(":", 1)[1].strip() for l in f if l.startswith("MemTotal")), None)
    except OSError:
        pass
    return dict(extended_precision_report=prec if prec is not None else xp.platform_precision(),
                float64=dict(eps=float(f64.eps), nmant=int(f64.nmant), mantissa_bits=int(f64.nmant) + 1),
                longdouble=dict(eps=float(fl.eps), nmant=int(fl.nmant), mantissa_bits=int(fl.nmant) + 1, itemsize=int(np.dtype(np.longdouble).itemsize), dtype=str(np.dtype(np.longdouble))),
                complex_longdouble_itemsize=int(np.dtype(np.clongdouble).itemsize), mpmath=dict(version=mpmath.__version__, dps=xp.MPMATH_DPS),
                machine=platform.machine(), platform=platform.platform(), cpu_model=cpu, cpu_count=os.cpu_count(), mem_total=mem,
                note="Tier 3 uses np.longdouble residual arithmetic only if >= 64 mantissa bits and the arithmetic probe passes; otherwise NOT EVALUABLE (never float64)")


def write_metadata(prob, out, R, comp, summ, cache_dir):
    import mpmath
    import scipy
    env = environment_info()
    ids = sorted({r["config"] for r in comp})
    meta = dict(
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"), command=" ".join([os.path.basename(sys.executable)] + sys.argv),
        problem=prob.name, source_commit=env.get("git_commit"), working_tree_dirty=env.get("git_dirty"),
        python=env.get("python"), numpy=np.__version__, scipy=scipy.__version__, mpmath=mpmath.__version__, platform=env.get("platform"),
        machine=env.get("machine"), cpu_count=env.get("cpu_count"), openblas_num_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
        precision=platform_precision_report(R.get("platform")),
        protocol=dict(file="docs/rgf_b3y_protocol.json", sha256=sha_file(PROTOCOL_PATH), commit=git_out("log", "-1", "--format=%H", "--", "docs/rgf_b3y_protocol.json"),
                      report_commit=git_out("log", "-1", "--format=%H", "--", "docs/RGF_B3Y_VALIDATION.md"), revision=PROTOCOL["revision_history"][-1] if PROTOCOL.get("revision_history") else None),
        base_commit_of_audit=PROTOCOL.get("base_commit"),
        potential_hashes={f"{x['config']}|1/{x['inv_dx']}|buffer{x['buffer_requested']}": x["potential_hash"] for x in R["identity"] if "potential_hash" in x},
        exact_energies=[dict(point=p["id"], config=p["config"], E=p["E"]) for p in prob.points.values()],
        discretization=dict(N_y=prob.Ny, N_q=prob.Nq, inv_dx_levels=list(prob.invs), buffer_requested=prob.buffer_req, buffer_control_requested=prob.buffer_ctrl_req),
        solver_settings=dict(M1="stationary.scattering_matrix (block recursion)", M2="splu permc_spec=COLAMD (default), Green blocks, smatrix_spec Fisher-Lee",
                             M3="splu permc_spec=MMD_AT_PLUS_A + 1 refinement step, S from the scattering state",
                             M4=dict(preconditioner_permc="MMD_AT_PLUS_A", refinement_max_iterations=xp.MAX_ITERATIONS, eta_max=xp.ETA_MAX, contraction_max=xp.CONTRACTION_MAX,
                                     fl_vs_state_tol=xp.FL_WF_TOL, ordering_tol=xp.ORDERING_TOL, error_indicator_max=xp.ERROR_INDICATOR_MAX, mpmath_dps=xp.MPMATH_DPS,
                                     chunk_bytes=xp.CHUNK_BYTES),
                             decomposition_M2="matrix sr.assemble_M, permc COLAMD, unit sources on open channels, smatrix_spec Fisher-Lee form",
                             decomposition_M3="matrix sw.assemble, permc MMD_AT_PLUS_A, method source, wave-function extraction"),
        limits=LIMITS, xp_gates=XPG,
        b3x_artifact_hashes=PROTOCOL.get("b3x_artifacts_sha256"), b3x_artifact_hashes_now={os.path.relpath(p, HERE): sha_file(p) for p in (B3X_CSV,)} if os.path.exists(B3X_CSV) else {},
        stopped_at=summ.get("stopped_at"), overall_status=summ.get("overall_status"), wall_time_s=summ.get("wall_time_s"),
        cache_dir_note="per-case caches live outside the result tree; the result files are assembled from them", **code_fingerprint("rgf_b3y"))
    with open(os.path.join(out, "run_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2, default=_json_default)


def record_tests(out, key="post_run_tests"):
    path = os.path.join(out, "run_metadata.json")
    with open(path) as f:
        m = json.load(f)
    runs = []
    for args in (["tests/test_rgf_b3y.py", "tests/test_rgf_b3y_driver.py"], ["tests"]):
        p = subprocess.run([sys.executable, "-m", "pytest", "-q"] + args, capture_output=True, text=True, cwd=HERE)
        last = [l for l in p.stdout.strip().splitlines() if l.strip()][-1]
        runs.append(dict(command="python -m pytest -q " + " ".join(args), exit_status=p.returncode, summary=last))
    m[key] = dict(recorded_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"), runs=runs,
                  note="executed after the artifacts were generated; written by validate_rgf_b3y.py --record-tests")
    with open(path, "w") as f:
        json.dump(m, f, indent=2)
    print(json.dumps(m[key], indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rgf_b3y_validation")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--record-tests", action="store_true")
    a = ap.parse_args()
    check_output_dir(a.out)
    if a.record_tests:
        record_tests(a.out)
        sys.exit(0)
    cache = a.cache or os.path.join(os.path.dirname(os.path.abspath(a.out)), "b3y_cache_" + ("syn" if a.synthetic else "real"))
    prob = SyntheticProblem() if a.synthetic else RealProblem()
    s = run_campaign(prob, a.out, cache, log=lambda m: print(m, flush=True), run_g0=not a.synthetic)
    print(json.dumps({k: s[k] for k in ("overall_status", "stopped_at", "wall_time_s")}, indent=1, default=_json_default))
