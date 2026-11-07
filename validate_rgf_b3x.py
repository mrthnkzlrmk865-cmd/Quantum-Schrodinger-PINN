r"""
B3X: extended independent validation of the RGF under fixed disorder
realizations. Protocol (fixed before any result): docs/RGF_B3X_VALIDATION.md
Part I and docs/rgf_b3x_protocol.json. Results: Part II of the same report.

Configurations: W in {100, 200} x zeta in {0.02, 0.05, 0.20}, seed 0, L_d = 1,
five energies each (kx = kx0 + j sigma_k, j = -2..2), dx = 1/512, N_y = 63, N_q = 127.
Three S-matrices are computed for every energy and compared as complex matrices:
  RGF                  stationary.scattering_matrix
  wave-function ref.   sparse_wavefunction (independent lead roots, potential blocks,
                       linear algebra and S assembly)
  Green-block ref.     sparse_reference (the B3 reference: independent solver, shared
                       lead roots and Fisher-Lee)
plus unitarity (both products), reciprocity (full, blocks, probabilities), open/closed
channel behaviour, and the refinements dx, N_q, buffer and N_y (at fixed N_q = 254).
A dx ladder 1/256 ... 1/2048 is computed as a DIAGNOSTIC only.

Usage: python validate_rgf_b3x.py [--out DIR] [--configs ids]
       python validate_rgf_b3x.py --record-tests [--out DIR]
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

import sparse_reference as sr
import sparse_wavefunction as sw
import validate_rgf_b3 as b3
from config import waveguide_config
from disorder import disorder_from_config
from potential import Potential
from provenance import code_fingerprint, environment_info
from stationary import build_device, scattering_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
PROTOCOL_PATH = os.path.join(HERE, "docs", "rgf_b3x_protocol.json")
with open(PROTOCOL_PATH) as _f:
    PROTOCOL = json.load(_f)
LIMITS = PROTOCOL["limits"]
SEED = PROTOCOL["configurations"]["seed"]
CONFIGS = {k: dict(W=v[0], zeta=v[1]) for k, v in PROTOCOL["configurations"]["ids"].items()}
BASE = dict(dx=PROTOCOL["base_discretization"]["dx"], Ny=PROTOCOL["base_discretization"]["N_y"],
            Nq=PROTOCOL["base_discretization"]["N_q"], buffer=PROTOCOL["base_discretization"]["buffer_requested"])
REF = PROTOCOL["refinements"]
VARIANTS = {
    "dx": dict(BASE, dx=REF["dx"]["variant"]),
    "Nq": dict(BASE, Nq=REF["N_q"]["variant"]),
    "buffer": dict(BASE, buffer=REF["buffer"]["variant"]),
    "Ny127_Nq254": dict(BASE, Ny=REF["N_y_fixed_Nq"]["variant"], Nq=REF["N_y_fixed_Nq"]["N_q"]),
    "Ny255": dict(BASE, Ny=255, Nq=511),
}
LADDER = tuple(REF["dx_ladder_diagnostic"])
assert BASE["dx"] == 1 / 512 and BASE["Nq"] == 127 and 2 * BASE["Nq"] == REF["N_q"]["variant"]
assert b3.BASE == BASE and b3.LIMITS["S_agreement"] == LIMITS["S_agreement"] and b3.LIMITS["dT1"] == LIMITS["dT1_refinement"]
assert b3.LIMITS["U"] == LIMITS["unitarity"] and b3.LIMITS["Q"] == LIMITS["reciprocity"]
B3_ID_MAP = {"A": "W100_z0.02", "B": "W100_z0.05", "C": "W100_z0.2", "D": "W200_z0.05"}   # existing B3 cases
passes = b3.passes


# ----------------------------------------------------------------------
def config_cfg(cid):
    p = CONFIGS[cid]
    return waveguide_config(disorder_strength=p["W"], disorder_corr_length=p["zeta"], slab_length=PROTOCOL["configurations"]["Ld"])


def config_potential(cid):
    cfg = config_cfg(cid)
    return cfg, Potential(cfg, disorder_from_config(cfg, seed=SEED))


def direct_samples(cid, x_start, dx, n_slices, Ly, Nq):
    """The same realization through a freshly built Potential/disorder object and own coordinates."""
    _, pot = config_potential(cid)
    xs = x_start + dx * np.arange(n_slices)
    y = np.arange(1, Nq + 1) * Ly / (Nq + 1)
    Xg, Yg = np.meshgrid(xs, y, indexing="ij")
    return pot.numpy(Xg, Yg)


def build_rgf(cid, cfg, pot, var):
    geo = b3.geometry(cfg, var["dx"], var["buffer"])
    rec = b3.Recorder(pot.numpy)
    dev = build_device(rec, geo["x_start"], geo["x_end"], var["dx"], cfg.Ly, var["Ny"], var["Nq"])
    ref = direct_samples(cid, geo["x_start"], var["dx"], dev.n_slices, cfg.Ly, var["Nq"])
    cmp = b3.compare_arrays(rec.calls[0], ref)
    return dev, geo, ref, dict(shape=list(ref.shape), potential_hash=b3.sha(ref), V_min=float(ref.min()), V_max=float(ref.max()),
                               V_lead_max=dev.V_lead_max, **cmp)


def offdiag_max(P):
    return float(np.max(P - np.diag(np.diag(P)))) if P.shape[0] > 1 else 0.0


def method_metrics(S, n_open):
    p = b3.sm.split_S(S, n_open)
    u, q = b3.unitarity(S), b3.reciprocity(S, n_open)
    return dict(U=u["U"], U_StS=u["U_StS"], U_SSt=u["U_SSt"], **{k: q[k] for k in q})


def compare_S(Sa, Sb, n_open):
    d = b3.s_difference(Sa, Sb)
    ia = np.unravel_index(int(np.argmax(d)), d.shape)
    Pa, Pb = np.abs(Sa) ** 2, np.abs(Sb) ** 2
    pa, pb = b3.sm.split_S(Pa, n_open), b3.sm.split_S(Pb, n_open)
    blocks = {k: float(np.max(np.abs(b3.sm.split_S(Sa, n_open)[k] - b3.sm.split_S(Sb, n_open)[k]))) for k in pa}
    off = max(float(np.max(np.abs((pa[k] - pb[k]) - np.diag(np.diag(pa[k] - pb[k]))))) if n_open > 1 else 0.0 for k in ("t_RL", "t_LR"))
    return dict(dS=float(d.max()), argmax_row=int(ia[0]), argmax_col=int(ia[1]), dS_rL=blocks["r_L"], dS_tLR=blocks["t_LR"],
                dS_tRL=blocks["t_RL"], dS_rR=blocks["r_R"], dP_all=float(np.max(np.abs(Pa - Pb))),
                dP_refl=float(max(np.max(np.abs(pa["r_L"] - pb["r_L"])), np.max(np.abs(pa["r_R"] - pb["r_R"])))),
                dP_trans=float(max(np.max(np.abs(pa["t_RL"] - pb["t_RL"])), np.max(np.abs(pa["t_LR"] - pb["t_LR"])))),
                dP_offdiag_trans=off, dT1=abs(b3.T1(Sa, n_open) - b3.T1(Sb, n_open)),
                dR1=abs(float(np.sum(np.abs(b3.sm.split_S(Sa, n_open)["r_L"][:, 0]) ** 2))
                        - float(np.sum(np.abs(b3.sm.split_S(Sb, n_open)["r_L"][:, 0]) ** 2))))


def closed_channel_checks(E, rgf_res, wf, Ep, dx=BASE["dx"]):
    """Open/closed channel behaviour (criteria C3, C4) and the closed-channel wave-function content."""
    lam_r, open_r, gam_r = rgf_res["modes"]["lam"], rgf_res["modes"]["open"], rgf_res["modes"]["gamma"]
    cl = ~open_r
    o = np.where(open_r)[0]
    eta = sw.retarded_root_eta(E, Ep, dx)
    # wave function of the RGF for the same sources, from its Green blocks (closed-channel amplitudes)
    G = rgf_res["G"]
    src = wf["source"][None, :]
    c1 = np.hstack([G["G11"][:, o] * src, G["G1N"][:, o] * src])
    cN = np.hstack([G["GN1"][:, o] * src, G["GNN"][:, o] * src])
    cl_idx = np.where(cl)[0]
    scale = float(max(np.max(np.abs(wf["c1"][cl_idx])), 1e-300)) if cl_idx.size else 0.0
    cdiff = float(max(np.max(np.abs(c1[cl_idx] - wf["c1"][cl_idx])), np.max(np.abs(cN[cl_idx] - wf["cN"][cl_idx])))) if cl_idx.size else 0.0
    return dict(
        n_closed=int(cl.sum()), gamma_closed_max_rgf=float(np.max(np.abs(gam_r[cl]))) if cl.any() else 0.0,
        gamma_closed_max_wf=float(np.max(np.abs(wf["gamma"][~wf["open"]]))) if (~wf["open"]).any() else 0.0,
        lam_closed_max_abs_rgf=float(np.max(np.abs(lam_r[cl]))), lam_closed_max_imag_rgf=float(np.max(np.abs(lam_r[cl].imag))),
        lam_closed_max_abs_wf=float(np.max(np.abs(wf["lam"][~wf["open"]]))), lam_closed_max_imag_wf=float(np.max(np.abs(wf["lam"][~wf["open"]].imag))),
        lam_open_max_dev_from_unit_circle=float(np.max(np.abs(np.abs(lam_r[open_r]) - 1.0))),
        lam_diff_rgf_vs_wf=float(np.max(np.abs(lam_r - wf["lam"]))), lam_diff_wf_vs_eta_root=float(np.max(np.abs(wf["lam"] - eta))),
        open_masks_equal=bool(np.array_equal(open_r, wf["open"])),
        closed_amp_diff_rgf_vs_wf=cdiff, closed_amp_scale=scale, closed_to_open_weight_slice1=wf["closed_to_open_weight_slice1"])


# ----------------------------------------------------------------------
def run_config(cid, log=print):
    cfg, pot = config_potential(cid)
    p = CONFIGS[cid]
    E_list = b3.case_energies(cfg)
    R = dict(energy=[], identity=[], refine=[], ladder=[], geometry={})
    dev, geo, ref, info = build_rgf(cid, cfg, pot, BASE)
    R["geometry"]["base"] = dict(BASE, **geo, V_lead_max=info["V_lead_max"], potential_hash=info["potential_hash"])

    def ident(vname, var, evaluator, info, geo):
        R["identity"].append(dict(config=cid, seed=SEED, W=p["W"], zeta=p["zeta"], variant=vname, evaluator=evaluator, dx=var["dx"],
                                  N_y=var["Ny"], N_q=var["Nq"], buffer=var["buffer"], n_slices=geo["n_slices"], x_start=geo["x_start"],
                                  x_end=geo["x_end"], disorder_region=list(cfg.disorder_region), n_disorder_modes=cfg.disorder_n_modes, **info))

    ident("base", BASE, "RGF device (build_device) vs direct Potential call", info, geo)
    Ep = dev.basis.E_perp
    base_T1, base_open = {}, {}
    t0 = time.time()
    for i, e in enumerate(E_list):
        E = e["E"]
        rg = scattering_matrix(E, dev)
        wf = sw.wavefunction_scattering_matrix(E, pot.numpy, geo["x_start"], BASE["dx"], dev.n_slices, cfg.Ly, BASE["Ny"], BASE["Nq"])
        gb = sr.sparse_scattering_matrix(E, pot.numpy, geo["x_start"], BASE["dx"], dev.n_slices, cfg.Ly, BASE["Ny"], BASE["Nq"])
        if i == 0:
            for lab, arr in (("wave-function reference (own sampling)", wf["V_samples"]), ("Green-block reference (own sampling)", gb["V_samples"])):
                ident("base", BASE, lab + " vs direct Potential call", dict(shape=list(arr.shape), potential_hash=b3.sha(arr), V_min=float(arr.min()),
                                                                           V_max=float(arr.max()), V_lead_max=info["V_lead_max"], **b3.compare_arrays(arr, ref)), geo)
        n_open = rg["n_open"]
        analytic = int(np.sum(0.5 * (np.arange(1, 5000) * np.pi / cfg.Ly) ** 2 < E))
        mm = {m: method_metrics(S, n_open) for m, S in (("rgf", rg["S"]), ("wf", wf["S"]), ("gb", gb["S"]))}
        row = dict(config=cid, seed=SEED, W=p["W"], zeta=p["zeta"], Ld=1.0, dx=BASE["dx"], N_y=BASE["Ny"], N_q=BASE["Nq"],
                   n_slices=dev.n_slices, matrix_order=dev.n_slices * BASE["Ny"], kx=e["kx"], E=E, n_open_rgf=n_open, n_open_wf=wf["n_open"],
                   n_open_gb=gb["n_open"], n_open_analytic=analytic, n_open_threshold_table=e["n_open_exact"], S_dim=int(rg["S"].shape[0]),
                   T1_rgf=b3.T1(rg["S"], n_open), R1_rgf=float(np.sum(np.abs(b3.sm.split_S(rg["S"], n_open)["r_L"][:, 0]) ** 2)),
                   offdiag_T_max_rgf=max(offdiag_max(np.abs(b3.sm.split_S(rg["S"], n_open)["t_RL"]) ** 2),
                                         offdiag_max(np.abs(b3.sm.split_S(rg["S"], n_open)["t_LR"]) ** 2)),
                   rgf_min_rcond=rg["min_rcond"], wf_residual=wf["residual_after_refinement"], wf_residual_before_refinement=wf["residual_before_refinement"],
                   gb_residual=gb["relative_residual"], wf_lu_nnz=wf["lu_nnz"])
        for k, v in compare_S(rg["S"], wf["S"], n_open).items():
            row[f"rgf_wf_{k}"] = v
        for k, v in compare_S(rg["S"], gb["S"], n_open).items():
            row[f"rgf_gb_{k}"] = v
        row["wf_gb_dS"] = float(np.max(b3.s_difference(wf["S"], gb["S"])))
        for m in ("rgf", "wf", "gb"):
            for k, v in mm[m].items():
                row[f"{k}_{m}"] = v
        row.update(closed_channel_checks(E, rg, wf, Ep))
        row["V_block_rel_diff_cosine_vs_rgf"] = float(np.max(np.abs(wf["V_blocks"] - dev.V)) / np.max(np.abs(wf["V_samples"])))
        row["V_block_rel_diff_cosine_vs_gb"] = float(np.max(np.abs(wf["V_blocks"] - gb["V_blocks"])) / np.max(np.abs(wf["V_samples"])))
        R["energy"].append(row)
        base_T1[i], base_open[i] = row["T1_rgf"], n_open
    log(f"  {cid} base (RGF + 2 references): {time.time() - t0:.0f}s, N = {dev.n_slices}")

    def run_variant(vname, var, role, name, ref_name, ref_vals):
        d2, g2, r2, inf2 = build_rgf(cid, cfg, pot, var)
        ident(vname, var, "RGF device (build_device) vs direct Potential call", inf2, g2)
        R["geometry"][vname] = dict(var, **g2, V_lead_max=inf2["V_lead_max"], potential_hash=inf2["potential_hash"])
        out, t0 = [], time.time()
        for i, e in enumerate(E_list):
            r = scattering_matrix(e["E"], d2)
            T = b3.T1(r["S"], r["n_open"])
            out.append(T)
            d = T - ref_vals[i]
            R["refine"].append(dict(config=cid, W=p["W"], zeta=p["zeta"], kx=e["kx"], E=e["E"], role=role, name=name, reference=ref_name,
                                    variant=vname, dx_ref=BASE["dx"], N_y_variant=var["Ny"], N_q_variant=var["Nq"], buffer_variant=g2["buffer_left"],
                                    n_slices_variant=g2["n_slices"], T1_reference=ref_vals[i], T1_variant=T, delta_T1=d, abs_delta_T1=abs(d),
                                    limit=LIMITS["dT1_refinement"], passed=passes(abs(d), LIMITS["dT1_refinement"]), n_open_variant=r["n_open"],
                                    U_variant=b3.unitarity(r["S"])["U"], V_lead_max=inf2["V_lead_max"]))
        log(f"  {cid} variant {vname}: N = {d2.n_slices}, N_y = {var['Ny']}: {time.time() - t0:.0f}s")
        return out

    T_dx = run_variant("dx", VARIANTS["dx"], "criterion", "dx", "base", base_T1)
    run_variant("Nq", VARIANTS["Nq"], "criterion", "Nq", "base", base_T1)
    run_variant("buffer", VARIANTS["buffer"], "criterion", "buffer", "base", base_T1)
    T_nq = [r["T1_variant"] for r in R["refine"] if r["variant"] == "Nq"]
    T_ny = run_variant("Ny127_Nq254", VARIANTS["Ny127_Nq254"], "criterion", "Ny_fixed_Nq254", "Nq254 (N_y = 63)", T_nq)
    if p["zeta"] == 0.02 and any(not r["passed"] for r in R["refine"] if r["name"] == "Ny_fixed_Nq254"):
        log(f"  {cid}: N_y criterion failed at 127 -> running N_y = 255 (protocol rule)")
        run_variant("Ny255", VARIANTS["Ny255"], "conditional_extension", "Ny_127_to_255", "Ny127_Nq254", T_ny)
    # diagnostic dx ladder (RGF only; cannot rescue a failed criterion)
    lad = {BASE["dx"]: [base_T1[i] for i in range(5)], VARIANTS["dx"]["dx"]: T_dx}
    for dxv in LADDER:
        if dxv in lad:
            continue
        d3, g3, r3, i3 = build_rgf(cid, cfg, pot, dict(BASE, dx=dxv))
        ident(f"ladder_dx_1_{round(1 / dxv)}", dict(BASE, dx=dxv), "RGF device (build_device) vs direct Potential call", i3, g3)
        t0 = time.time()
        lad[dxv] = [b3.T1(*(lambda r: (r["S"], r["n_open"]))(scattering_matrix(e["E"], d3))) for e in E_list]
        log(f"  {cid} ladder dx = 1/{round(1 / dxv)}: N = {d3.n_slices}: {time.time() - t0:.0f}s")
    for i, e in enumerate(E_list):
        t = [lad[dxv][i] for dxv in LADDER]
        o1 = float(np.log2(abs(t[0] - t[1]) / abs(t[1] - t[2]))) if t[0] != t[1] and t[1] != t[2] else float("nan")
        o2 = float(np.log2(abs(t[1] - t[2]) / abs(t[2] - t[3]))) if t[1] != t[2] and t[2] != t[3] else float("nan")
        ok = bool(1.8 <= o1 <= 2.2 and 1.8 <= o2 <= 2.2)
        ext = t[3] + (t[3] - t[2]) / 3.0 if ok else float("nan")
        R["ladder"].append(dict(config=cid, W=p["W"], zeta=p["zeta"], kx=e["kx"], E=e["E"], label="diagnostic",
                                **{f"T1_dx_1_{round(1 / d)}": v for d, v in zip(LADDER, t)},
                                delta_512_to_1024=t[2] - t[1], delta_1024_to_2048=t[3] - t[2], order_256_512_1024=o1, order_512_1024_2048=o2,
                                richardson_applied=ok, T1_richardson=ext, est_error_dx512=(t[1] - ext) if ok else float("nan"),
                                est_error_dx1024=(t[2] - ext) if ok else float("nan")))
    return R


# ----------------------------------------------------------------------
def _crit(rows, key, limit):
    return b3._crit([dict(r, case=r["config"]) for r in rows], key, limit)


def summarize(R):
    en = R["energy"]
    s = {}
    for r in en:
        r["_tr"] = max(r["rgf_wf_dP_all"], r["rgf_gb_dP_all"], r["rgf_wf_dT1"], r["rgf_wf_dR1"], r["rgf_gb_dT1"], r["rgf_gb_dR1"])
        r["_off"] = max(r["rgf_wf_dP_offdiag_trans"], r["rgf_gb_dP_offdiag_trans"])
        r["_U"] = max(r["U_rgf"], r["U_wf"], r["U_gb"])
        r["_Q"] = max(r["Q_rgf"], r["Q_wf"], r["Q_gb"])
    ident = R["identity"]
    s["C1_same_realization_bitwise"] = dict(
        passed=bool(all(r["bitwise_equal"] and r["coords_bitwise_equal"] if "coords_bitwise_equal" in r else r["bitwise_equal"] for r in ident)),
        n_checks=len(ident), max_abs_diff=float(max(r["max_abs_diff"] for r in ident)),
        max_V_lead=float(max(r["V_lead_max"] for r in ident)))
    v = float(max(max(r["V_block_rel_diff_cosine_vs_rgf"], r["V_block_rel_diff_cosine_vs_gb"]) for r in en))
    s["C2_potential_blocks_independent_routes"] = dict(value=v, limit=LIMITS["potential_block_agreement_relative"], passed=passes(v, LIMITS["potential_block_agreement_relative"]))
    mism = sum(not (r["n_open_rgf"] == r["n_open_wf"] == r["n_open_gb"] == r["n_open_analytic"] == r["n_open_threshold_table"]) for r in en)
    s["C3_open_channel_count"] = dict(value=int(mism), limit=0, passed=mism == 0, n_open_range=[int(min(r["n_open_rgf"] for r in en)), int(max(r["n_open_rgf"] for r in en))])
    c4 = dict(max_gamma_closed=float(max(max(r["gamma_closed_max_rgf"], r["gamma_closed_max_wf"]) for r in en)),
              max_closed_lambda_imag=float(max(max(r["lam_closed_max_imag_rgf"], r["lam_closed_max_imag_wf"]) for r in en)),
              max_closed_lambda_abs=float(max(max(r["lam_closed_max_abs_rgf"], r["lam_closed_max_abs_wf"]) for r in en)),
              min_closed_channels=int(min(r["n_closed"] for r in en)), max_closed_channels=int(max(r["n_closed"] for r in en)),
              max_lambda_diff_rgf_vs_wf=float(max(r["lam_diff_rgf_vs_wf"] for r in en)),
              max_lambda_diff_wf_vs_eta_root=float(max(r["lam_diff_wf_vs_eta_root"] for r in en)),
              open_masks_equal=bool(all(r["open_masks_equal"] for r in en)),
              report_only_closed_channel_amplitude_diff_rgf_vs_wf=float(max(r["closed_amp_diff_rgf_vs_wf"] for r in en)),
              report_only_closed_channel_amplitude_scale=float(max(r["closed_amp_scale"] for r in en)),
              report_only_closed_to_open_weight_slice1_max=float(max(r["closed_to_open_weight_slice1"] for r in en)))
    c4["passed"] = bool(c4["max_gamma_closed"] == 0.0 and c4["max_closed_lambda_imag"] == 0.0 and c4["max_closed_lambda_abs"] < 1.0
                        and c4["min_closed_channels"] >= 1 and c4["open_masks_equal"]
                        and passes(max(c4["max_lambda_diff_rgf_vs_wf"], c4["max_lambda_diff_wf_vs_eta_root"]), LIMITS["retarded_root_agreement"]))
    s["C4_closed_channels"] = c4
    for key, col, lim in (("C5_S_agreement_RGF_vs_wavefunction", "rgf_wf_dS", LIMITS["S_agreement"]),
                          ("C6_S_agreement_RGF_vs_greenblock", "rgf_gb_dS", LIMITS["S_agreement"]),
                          ("C7_T_R_agreement", "_tr", LIMITS["T_R_agreement"]),
                          ("C8_offdiag_transmission_agreement", "_off", LIMITS["offdiag_transmission_agreement"]),
                          ("C9_unitarity_all_methods", "_U", LIMITS["unitarity"]),
                          ("C10_reciprocity_all_methods", "_Q", LIMITS["reciprocity"])):
        s[key] = _crit(en, col, lim)
    s["C5_S_agreement_RGF_vs_wavefunction"]["S_dims"] = sorted({r["S_dim"] for r in en})
    s["C5_S_agreement_RGF_vs_wavefunction"]["secondary_wf_vs_greenblock_max"] = float(max(r["wf_gb_dS"] for r in en))
    s["C7_T_R_agreement"]["parts"] = {k: float(max(r[k] for r in en)) for k in ("rgf_wf_dP_all", "rgf_wf_dP_refl", "rgf_wf_dP_trans", "rgf_wf_dT1", "rgf_wf_dR1",
                                                                              "rgf_gb_dP_all", "rgf_gb_dT1", "rgf_gb_dR1")}
    s["C9_unitarity_all_methods"]["parts"] = {m: dict(StS=float(max(r[f"U_StS_{m}"] for r in en)), SSt=float(max(r[f"U_SSt_{m}"] for r in en))) for m in ("rgf", "wf", "gb")}
    s["C10_reciprocity_all_methods"]["parts"] = {m: {k: float(max(r[f"{k}_{m}"] for r in en)) for k in ("Q_full", "Q_rL", "Q_rR", "Q_t", "Q_prob")} for m in ("rgf", "wf", "gb")}
    s["C10_reciprocity_all_methods"]["report_only_same_direction_asymmetry_max"] = float(max(max(r[f"same_dir_RL_report_only_{m}"] for m in ("rgf", "wf", "gb")) for r in en))
    s["report_only_channel_mixing_max_offdiag_T_rgf"] = float(max(r["offdiag_T_max_rgf"] for r in en))
    ref = R["refine"]
    for key, name in (("C11_dx_refinement", "dx"), ("C12_Nq_refinement", "Nq"), ("C13_buffer_refinement", "buffer"), ("C14_Ny_refinement_fixed_Nq254", "Ny_fixed_Nq254")):
        rows = [r for r in ref if r["name"] == name and r["role"] == "criterion"]
        s[key] = _crit(rows, "abs_delta_T1", LIMITS["dT1_refinement"])
        s[key]["rows_failed"] = [dict(config=r["config"], kx=r["kx"], E=r["E"], abs_delta_T1=r["abs_delta_T1"]) for r in rows if not r["passed"]]
        s[key]["n_rows"] = len(rows)
    ext = [r for r in ref if r["role"] == "conditional_extension"]
    s["conditional_Ny255"] = dict(triggered=bool(ext), rule=PROTOCOL["refinements"]["N_y_255"],
                                  rows=[dict(config=r["config"], E=r["E"], abs_delta_T1=r["abs_delta_T1"], passed=r["passed"]) for r in ext],
                                  max_abs_delta_127_to_255=float(max([r["abs_delta_T1"] for r in ext], default=float("nan"))),
                                  status="labelled extension; does not replace the C14 result")
    lad = R["ladder"]
    ok = [r for r in lad if r["richardson_applied"]]
    s["dx_ladder_diagnostic"] = dict(
        label="DIAGNOSTIC (not a criterion)", n=len(lad), n_richardson_applied=len(ok),
        orders=dict(min=float(np.nanmin([r[k] for r in lad for k in ("order_256_512_1024", "order_512_1024_2048")])),
                    max=float(np.nanmax([r[k] for r in lad for k in ("order_256_512_1024", "order_512_1024_2048")]))),
        max_abs_est_error_dx512_per_config={c: float(np.nanmax([abs(r["est_error_dx512"]) for r in lad if r["config"] == c] or [float("nan")])) for c in sorted({r["config"] for r in lad})},
        max_abs_est_error_dx1024_per_config={c: float(np.nanmax([abs(r["est_error_dx1024"]) for r in lad if r["config"] == c] or [float("nan")])) for c in sorted({r["config"] for r in lad})},
        max_abs_delta_1024_to_2048_per_config={c: float(max(abs(r["delta_1024_to_2048"]) for r in lad if r["config"] == c)) for c in sorted({r["config"] for r in lad})})
    s["criteria_passed"] = {k: bool(v["passed"]) for k, v in s.items() if isinstance(v, dict) and k.startswith("C") and "passed" in v}
    assert len(s["criteria_passed"]) == 14, s["criteria_passed"].keys()
    s["failed_criteria"] = [k for k, v in s["criteria_passed"].items() if not v]
    s["all_criteria_evaluated"] = True
    s["b3x_full_pass"] = not s["failed_criteria"]
    s["overall_B3X_status"] = "PASS" if s["b3x_full_pass"] else "FAIL: " + ", ".join(s["failed_criteria"])
    return s


def per_config_summary(R, s):
    rows = []
    for cid, p in CONFIGS.items():
        en = [r for r in R["energy"] if r["config"] == cid]
        if not en:
            continue
        ref = [r for r in R["refine"] if r["config"] == cid and r["role"] == "criterion"]
        lad = [r for r in R["ladder"] if r["config"] == cid]
        row = dict(config=cid, W=p["W"], zeta=p["zeta"], seed=SEED, n_energies=len(en), n_open_range=f"{min(r['n_open_rgf'] for r in en)}-{max(r['n_open_rgf'] for r in en)}",
                   closed_channels_range=f"{min(r['n_closed'] for r in en)}-{max(r['n_closed'] for r in en)}",
                   max_dS_rgf_wf=max(r["rgf_wf_dS"] for r in en), max_dS_rgf_gb=max(r["rgf_gb_dS"] for r in en), max_dS_wf_gb=max(r["wf_gb_dS"] for r in en),
                   max_TR_discrepancy=max(r["_tr"] for r in en), max_offdiag_trans_discrepancy=max(r["_off"] for r in en),
                   max_U=max(r["_U"] for r in en), max_Q=max(r["_Q"] for r in en), max_offdiag_T_rgf=max(r["offdiag_T_max_rgf"] for r in en),
                   min_rgf_rcond=min(r["rgf_min_rcond"] for r in en))
        for nm, k in (("dx", "dx"), ("Nq", "Nq"), ("buffer", "buffer"), ("Ny_fixed_Nq254", "Ny_fixed_Nq254")):
            row[f"max_abs_dT1_{nm}"] = max(r["abs_delta_T1"] for r in ref if r["name"] == k)
        row["diag_max_abs_est_error_dx512"] = float(np.nanmax([abs(r["est_error_dx512"]) for r in lad] or [float("nan")]))
        row["diag_max_abs_est_error_dx1024"] = float(np.nanmax([abs(r["est_error_dx1024"]) for r in lad] or [float("nan")]))
        for key, col, lim in (("S_rgf_wf", "max_dS_rgf_wf", LIMITS["S_agreement"]), ("S_rgf_gb", "max_dS_rgf_gb", LIMITS["S_agreement"]),
                              ("TR", "max_TR_discrepancy", LIMITS["T_R_agreement"]), ("offdiag", "max_offdiag_trans_discrepancy", LIMITS["offdiag_transmission_agreement"]),
                              ("unitarity", "max_U", LIMITS["unitarity"]), ("reciprocity", "max_Q", LIMITS["reciprocity"]),
                              ("dx", "max_abs_dT1_dx", LIMITS["dT1_refinement"]), ("Nq", "max_abs_dT1_Nq", LIMITS["dT1_refinement"]),
                              ("buffer", "max_abs_dT1_buffer", LIMITS["dT1_refinement"]), ("Ny", "max_abs_dT1_Ny_fixed_Nq254", LIMITS["dT1_refinement"])):
            row[f"pass_{key}"] = passes(row[col], lim)
        rows.append(row)
    return rows


def cross_check_b3(R):
    path = os.path.join(HERE, "results", "rgf_b3_validation", "rgf_sparse_comparison.csv")
    if not os.path.exists(path):
        return []
    old = list(csv.DictReader(open(path)))
    rows = []
    for o in old:
        cid = B3_ID_MAP[o["case"]]
        for r in R["energy"]:
            if r["config"] == cid and r["E"] == float(o["E"]):
                rows.append(dict(config=cid, b3_case=o["case"], E=r["E"], T1_b3x_rgf=r["T1_rgf"], T1_b3_rgf=float(o["T1_RGF"]), T1_b3_sparse=float(o["T1_sparse"]),
                                 abs_diff_rgf=abs(r["T1_rgf"] - float(o["T1_RGF"])), abs_diff_b3x_rgf_vs_b3_sparse=abs(r["T1_rgf"] - float(o["T1_sparse"])),
                                 label="report-only cross-check with the existing B3"))
    return rows


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r if not k.startswith("_")))
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)


def sha_file(path):
    with open(path, "rb") as f:
        return "sha256:" + hashlib.sha256(f.read()).hexdigest()


def git_out(*args):
    try:
        return subprocess.run(["git", *args], capture_output=True, text=True, cwd=HERE, timeout=20).stdout.strip() or None
    except Exception:                                   # noqa: BLE001
        return None


def record_tests(out, key="post_run_tests"):
    path = os.path.join(out, "run_metadata.json")
    with open(path) as f:
        m = json.load(f)
    runs = []
    for args in (["tests/test_rgf_b3x.py"], ["tests"]):
        p = subprocess.run([sys.executable, "-m", "pytest", "-q"] + args, capture_output=True, text=True, cwd=HERE)
        last = [l for l in p.stdout.strip().splitlines() if l.strip()][-1]
        runs.append(dict(command="python -m pytest -q " + " ".join(args), exit_status=p.returncode, summary=last))
    m[key] = dict(recorded_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"), runs=runs,
                  note="executed after the artifacts were generated; written by validate_rgf_b3x.py --record-tests")
    with open(path, "w") as f:
        json.dump(m, f, indent=2)
    print(json.dumps(m[key], indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rgf_b3x_validation")
    ap.add_argument("--configs", default=",".join(CONFIGS))
    ap.add_argument("--record-tests", action="store_true")
    a = ap.parse_args()
    if a.record_tests:
        record_tests(a.out)
        sys.exit(0)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    R = dict(energy=[], identity=[], refine=[], ladder=[], geometry={})
    for cid in a.configs.split(","):
        print(f"config {cid}: {CONFIGS[cid]}", flush=True)
        r = run_config(cid, log=lambda m: print(m, flush=True))
        for k in ("energy", "identity", "refine", "ladder"):
            R[k] += r[k]
        R["geometry"][cid] = r["geometry"]
    s = summarize(R)
    pc = per_config_summary(R, s)
    xc = cross_check_b3(R)
    write_csv(os.path.join(a.out, "per_energy_comparison.csv"), R["energy"])
    write_csv(os.path.join(a.out, "per_configuration_summary.csv"), pc)
    write_csv(os.path.join(a.out, "refinement_convergence.csv"), R["refine"])
    write_csv(os.path.join(a.out, "dx_ladder_diagnostic.csv"), R["ladder"])
    write_csv(os.path.join(a.out, "potential_identity.csv"), R["identity"])
    if xc:
        write_csv(os.path.join(a.out, "cross_check_b3.csv"), xc)
        s["cross_check_with_existing_B3"] = dict(label="report-only", n=len(xc), max_abs_diff_rgf_T1=max(r["abs_diff_rgf"] for r in xc),
                                                 max_abs_diff_vs_b3_sparse_T1=max(r["abs_diff_b3x_rgf_vs_b3_sparse"] for r in xc))
    s.update(wall_time_s=time.time() - t0, limits=LIMITS, scope="numerical-method validation of fixed disorder realizations; not a statistical result",
             configurations=CONFIGS, seed=SEED, base_discretization=BASE)
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(s, f, indent=2)
    env = environment_info()
    import scipy
    from disorder import disorder_from_config as dfc
    meta = dict(
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"), command=" ".join([os.path.basename(sys.executable)] + sys.argv),
        source_commit=env.get("git_commit"), working_tree_dirty=env.get("git_dirty"), python=env.get("python"), numpy=np.__version__, scipy=scipy.__version__,
        platform=env.get("platform"), machine=env.get("machine"), cpu_count=env.get("cpu_count"), openblas_num_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
        protocol=dict(file="docs/rgf_b3x_protocol.json", sha256=sha_file(PROTOCOL_PATH), commit=git_out("log", "-1", "--format=%H", "--", "docs/rgf_b3x_protocol.json"),
                      report_commit=git_out("log", "-1", "--format=%H", "--", "docs/RGF_B3X_VALIDATION.md")),
        random_seed=SEED, configurations={c: dict(CONFIGS[c], seed=SEED, Ld=1.0, disorder_n_modes=config_cfg(c).disorder_n_modes, disorder_edge=config_cfg(c).disorder_edge,
                                                  disorder_region=list(config_cfg(c).disorder_region),
                                                  disorder_params_hash=b3.sha(np.concatenate([(d := dfc(config_cfg(c), seed=SEED)).qx, d.qy, d.phi])))
                                          for c in a.configs.split(",")},
        energies=[dict(config=r["config"], kx=r["kx"], E=r["E"], n_open=r["n_open_rgf"]) for r in R["energy"]],
        base_discretization=BASE, variants=VARIANTS, dx_ladder=list(LADDER), geometry=R["geometry"],
        methods=dict(rgf="stationary.scattering_matrix -> rgf_endpoint_blocks (block recursion) + smatrix_spec Fisher-Lee",
                     wavefunction_reference="sparse_wavefunction.wavefunction_scattering_matrix: own lead roots, own cosine-sum potential blocks, "
                                            "scipy.sparse.linalg.splu (SuperLU, MMD_AT_PLUS_A) + 1 refinement step, S from scattering state",
                     greenblock_reference="sparse_reference.sparse_scattering_matrix: splu (COLAMD), Green blocks + smatrix_spec Fisher-Lee"),
        limits=LIMITS, wall_time_s=s["wall_time_s"], **code_fingerprint("rgf_b3x"))
    with open(os.path.join(a.out, "run_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps({k: s[k] for k in ("criteria_passed", "failed_criteria", "overall_B3X_status", "wall_time_s")}, indent=1))
