r"""
EXPERIMENTAL diagnostic (single point): original RGF vs the small-quantity W/Z formulation at the worst B1 point
(configuration F, E = 125.69231487269683, channel n = 5). See PREDICTIONS.md (committed before this run) and REPORT.md.

Production code is only IMPORTED read-only (validate_rgf_b1.device / spec_energies, stationary.scattering_matrix) to obtain the
committed device, the production S-matrix and its double lead data. Nothing outside experiments/wz_rgf_diagnostic/ is written.

Exact references: mpmath at DPS digits, closed form (Chebyshev U_k of the Toeplitz block with corrected corners) for every
represented lattice, cross-checked by an independent mpmath recursion.
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import csv
import datetime as _dt
import hashlib
import json
import subprocess
import sys
import time

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from stationary import scattering_matrix          # noqa: E402  (read-only use)
from validate_rgf_b1 import device, spec_energies  # noqa: E402  (read-only use)
import wz_rgf as wz                                 # noqa: E402

DPS = 110
E_TARGET = "125.69231487269683"
CHANNEL = 5
OUT = os.path.join(HERE, "results")


# ----------------------------------------------------------------------
# exact scalar S for a represented lattice
# ----------------------------------------------------------------------
def mpc(z):
    return mp.mpc(complex(z).real, complex(z).imag)


def exact_scalar(a, t, sigma, gamma, N):
    """Scalar clean lattice: diagonal a on every slice, off-diagonal +t, self-energy sigma added (subtracted from M) at slices 1 and N,
    broadening gamma (both leads identical). Closed form: M = T + s (e_1 e_1^T + e_N e_N^T), s = -sigma,
    det M = D_N + 2 s D_{N-1} + s^2 D_{N-2}, D_k = t^k U_k(a / 2t);  G_N1 = (-1)^(N+1) t^(N-1) / det;  G_11 = (D_{N-1} + s D_{N-2}) / det."""
    a, t, s, gamma = mp.mpf(a), mp.mpf(t), -sigma, mp.mpf(gamma)
    x = a / (2 * t)
    rho = x + mp.sqrt(x * x - 1)
    U = lambda k: (rho ** (k + 1) - rho ** (-(k + 1))) / (rho - 1 / rho)
    mu = s / t
    det_n = U(N) + 2 * mu * U(N - 1) + mu * mu * U(N - 2)                      # det M / t^N
    GN1 = (-1) ** (N + 1) / (t * det_n)
    G11 = (U(N - 1) + mu * U(N - 2)) / (t * det_n)
    tt = 1j * gamma * GN1
    r = 1j * gamma * G11 - 1
    return dict(t_RL=tt, t_LR=tt, r_L=r, r_R=r, GN1=GN1, G11=G11)


def exact_scalar_recursion(a, t, sigma, gamma, N):
    """Independent cross-check: the original recursion carried out at DPS digits."""
    a, t = mp.mpf(a), mp.mpf(t)
    g = GN1 = G1N = G11 = None
    for j in range(N):
        D = mp.mpc(a)
        if j == 0:
            D -= sigma
        if j == N - 1:
            D -= sigma
        if j > 0:
            D -= t * t * g
        gn = 1 / D
        if j == 0:
            g = GN1 = G1N = G11 = gn
            continue
        G11 = G11 + t * t * G1N * gn * GN1
        G1N = -t * G1N * gn
        GN1 = -t * gn * GN1
        g = gn
    gamma = mp.mpf(gamma)
    return dict(t_RL=1j * gamma * GN1, t_LR=1j * gamma * G1N, r_L=1j * gamma * G11 - 1, r_R=1j * gamma * g - 1)


def root_open(c):
    """Retarded root of lambda + 1/lambda = 2c for |c| < 1 (open): lambda = c + i sqrt(1 - c^2)."""
    return c + 1j * mp.sqrt(1 - c * c)


def lattice(name, a, t, sigma, gamma, note):
    a, t = mp.mpf(a), mp.mpf(t)
    c_dev = -a / (2 * t)
    lam = -sigma / t
    c_lead = (lam + 1 / lam) / 2
    return dict(name=name, a=a, t=t, sigma=sigma, gamma=mp.mpf(gamma), c_dev=c_dev, lam=lam, c_lead=c_lead,
                c_mismatch=mp.mpf(abs(c_lead - c_dev)), gamma_vs_2tImlam=mp.mpf(abs(mp.mpf(gamma) - 2 * t * mp.im(lam))),
                consistent=bool(abs(c_lead - c_dev) < mp.mpf(10) ** (-DPS + 10) and abs(mp.mpf(gamma) - 2 * t * mp.im(lam)) < mp.mpf(10) ** (-DPS + 10)),
                note=note)


# ----------------------------------------------------------------------
def metrics(t_num, t_ref):
    t_num, t_ref = mpc(t_num) if not isinstance(t_num, mp.mpc) else t_num, t_ref
    q = (t_num - t_ref) / t_ref
    return dict(q_re=float(mp.re(q)), q_im=float(mp.im(q)), abs_q=float(abs(q)),
                eps_A=float(abs(t_num) / abs(t_ref) - 1), eps_phi=float(mp.arg(t_num * mp.conj(t_ref))))


def S2(d):
    return np.array([[complex(d["r_L"]), complex(d["t_LR"])], [complex(d["t_RL"]), complex(d["r_R"])]])


def S2mp(d):
    return mp.matrix([[d["r_L"], d["t_LR"]], [d["t_RL"], d["r_R"]]])


def unitarity_defect_mp(S):
    D = S.H * S - mp.eye(2)
    D2 = S * S.H - mp.eye(2)
    return float(max(max(abs(D[i, j]) for i in range(2) for j in range(2)), max(abs(D2[i, j]) for i in range(2) for j in range(2))))


def git(*a):
    try:
        return subprocess.run(["git", *a], capture_output=True, text=True, cwd=ROOT, timeout=20).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def main():
    mp.mp.dps = DPS
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    # ---- committed device and energy ------------------------------------------------------------
    dev = device(63, 1 / 1024, 3.0)
    Es = spec_energies(dev.basis.E_perp)
    E = [e for e in Es if repr(float(e)) == E_TARGET]
    assert len(E) == 1, "target energy not on the committed spec grid"
    E = float(E[0])
    n = CHANNEL - 1
    Ep = float(dev.basis.E_perp[n])
    dx, N, t = float(dev.dx), int(dev.n_slices), float(dev.t_x)
    assert N == 3483 and dx == 1 / 1024 and t == 2.0 ** 19
    # ---- production solver (unchanged code, read-only) -------------------------------------------
    r = scattering_matrix(E, dev)
    m = r["modes"]
    k = list(np.where(m["open"])[0]).index(n)
    prod = dict(t_RL=r["t_RL"][k, k], t_LR=r["t_LR"][k, k], r_L=r["r_L"][k, k], r_R=r["r_R"][k, k])
    lam_code, gam_code, sig_code = complex(m["lam"][n]), float(m["gamma"][n]), complex(m["sigma"][n])
    off_t = float(max(np.max(np.abs(r["t_RL"] - np.diag(np.diag(r["t_RL"])))), np.max(np.abs(r["r_L"] - np.diag(np.diag(r["r_L"]))))))
    # ---- stored inputs ------------------------------------------------------------------------------
    a_dev = E - (2.0 * t + Ep)                                   # fl(E - fl(2t + E_perp)), as D = E*I - H in stationary.py
    delta_fl = (E - Ep) * dx * dx                                # W/Z representation
    # ---- represented lattices ------------------------------------------------------------------------
    mE, mEp, mdx, mt = mp.mpf(E), mp.mpf(Ep), mp.mpf(dx), mp.mpf(t)
    c_L = 1 - (mE - mEp) * mdx ** 2
    lamL = root_open(c_L)
    c_dev = -mp.mpf(a_dev) / (2 * mt)
    lam_dev = root_open(c_dev)
    c_W = 1 - mp.mpf(delta_fl)
    lam_W = root_open(c_W)
    lat = {
        "L": lattice("L", -2 * mt * c_L, mt, -mt * lamL, 2 * mt * mp.im(lamL), "exact lattice; double E, E_perp taken as exact"),
        "R_dev": lattice("R_dev", mp.mpf(a_dev), mt, -mt * lam_dev, 2 * mt * mp.im(lam_dev), "stored device diagonal a_dev; exact leads consistent with c_dev"),
        "R_code": lattice("R_code", mp.mpf(a_dev), mt, mpc(sig_code), mp.mpf(gam_code), "stored device diagonal a_dev; production double lead data (sigma, Gamma from smatrix_spec.lead_modes)"),
        "R_W": lattice("R_W", mt * (2 * mp.mpf(delta_fl) - 2), mt, -mt * lam_W, 2 * mt * mp.im(lam_W), "W/Z device Delta_fl = fl(E - E_perp) dx^2; exact leads consistent with c_W"),
        "R_W_code_lead": lattice("R_W_code_lead", mt * (2 * mp.mpf(delta_fl) - 2), mt, mpc(sig_code), mp.mpf(gam_code), "W/Z device Delta_fl; production double lead data"),
    }
    ref, xcheck = {}, {}
    for key, L in lat.items():
        ref[key] = exact_scalar(L["a"], L["t"], L["sigma"], L["gamma"], N)
        rec = exact_scalar_recursion(L["a"], L["t"], L["sigma"], L["gamma"], N)
        xcheck[key] = float(max(abs(ref[key][q] - rec[q]) for q in ("t_RL", "t_LR", "r_L", "r_R")))
    # analytic identity for consistent leads: t = lambda^(N-1)
    analytic_check = {kk: float(abs(ref[kk]["t_RL"] - lam ** (N - 1))) for kk, lam in (("L", lamL), ("R_dev", lam_dev), ("R_W", lam_W))}
    # ---- double solvers ---------------------------------------------------------------------------------
    rep = wz.original_scalar_recurrence(a_dev, sig_code, t, N, record_d=True)
    replica = dict(t_RL=1j * gam_code * rep["GN1"], t_LR=1j * gam_code * rep["G1N"], r_L=1j * gam_code * rep["G11"] - 1, r_R=1j * gam_code * rep["GNN"] - 1)
    Delta = np.full((N, 1, 1), delta_fl)
    lead_d = wz.lead_from_delta(np.array([delta_fl]))
    lead_c = wz.lead_from_code_data(np.array([lam_code]), np.array([gam_code]), t)
    out_wz = {}
    for tag, ld in (("WZ_delta_lead", lead_d), ("WZ_code_lead", lead_c)):
        G = wz.wz_endpoint_blocks(Delta, ld["lam_minus_1"], record_w=True)
        S = wz.wz_fisher_lee(G, ld["sin_theta"], ld["open"])
        out_wz[tag] = dict(t_RL=S["t_RL"][0, 0], t_LR=S["t_LR"][0, 0], r_L=S["r_L"][0, 0], r_R=S["r_R"][0, 0],
                           W_max=G["W_max"], Z=complex(G["Z"][0, 0]), w=[complex(w[0, 0]) for w in G["w_trace"]])
    solvers = {"original_production": prod, "original_scalar_replica": replica,
               "WZ_delta_lead": out_wz["WZ_delta_lead"], "WZ_code_lead": out_wz["WZ_code_lead"]}
    # ---- comparisons ------------------------------------------------------------------------------------
    comps = [("C1", "original_production", "L", "total error of the original formulation"),
             ("C2", "original_production", "R_code", "recurrence rounding of the original (same represented lattice)"),
             ("C3", "WZ_delta_lead", "L", "total error of W/Z with delta-based leads"),
             ("C4", "WZ_delta_lead", "R_W", "W/Z rounding (recurrence + double delta-based lead data) on its represented lattice"),
             ("C5", "WZ_code_lead", "L", "total error of W/Z with the production lead data"),
             ("C6", "WZ_code_lead", "R_W_code_lead", "W/Z recurrence rounding on its represented lattice"),
             ("X1", "original_scalar_replica", "R_code", "scalar replica of the original recurrence (trajectory source)")]
    rows = []
    for cid, s_name, r_name, what in comps:
        sv, rf = solvers[s_name], ref[r_name]
        mt_ = metrics(mpc(sv["t_RL"]), rf["t_RL"])
        ml_ = metrics(mpc(sv["t_LR"]), rf["t_LR"])
        Sd = S2(sv)
        U = float(max(np.max(np.abs(Sd.conj().T @ Sd - np.eye(2))), np.max(np.abs(Sd @ Sd.conj().T - np.eye(2)))))
        rows.append(dict(id=cid, solver=s_name, reference=r_name, what=what,
                         eps_phi_tRL=mt_["eps_phi"], eps_phi_tLR=ml_["eps_phi"], eps_A_tRL=mt_["eps_A"], eps_A_tLR=ml_["eps_A"],
                         abs_q_tRL=mt_["abs_q"], q_re_tRL=mt_["q_re"], q_im_tRL=mt_["q_im"],
                         abs_r_solver=float(max(abs(complex(sv["r_L"])), abs(complex(sv["r_R"])))),
                         abs_r_reference=float(max(abs(rf["r_L"]), abs(rf["r_R"]))),
                         abs_dr_L=float(abs(mpc(sv["r_L"]) - rf["r_L"])),
                         unitarity_defect_solver=U, reciprocity_defect_solver=float(abs(complex(sv["t_LR"]) - complex(sv["t_RL"])))))
    # reference-to-reference differences (pure lattice effects)
    refdiff = []
    for a_, b_, what in (("R_dev", "L", "device representation (original)"), ("R_code", "R_dev", "lead mismatch (original)"),
                         ("R_W", "L", "device representation (W/Z)"), ("R_W_code_lead", "R_W", "lead mismatch (W/Z with production leads)")):
        mt_ = metrics(ref[a_]["t_RL"], ref[b_]["t_RL"])
        refdiff.append(dict(lattice=a_, minus=b_, what=what, eps_phi_tRL=mt_["eps_phi"], eps_A_tRL=mt_["eps_A"], abs_q_tRL=mt_["abs_q"],
                            abs_dr_L=float(abs(ref[a_]["r_L"] - ref[b_]["r_L"])), abs_r_lattice=float(abs(ref[a_]["r_L"])),
                            unitarity_defect_lattice=unitarity_defect_mp(S2mp(ref[a_]))))
    # ---- telescoping identities (complex t_RL and r_L) -----------------------------------------------------
    def tele(parts, total):
        return float(abs(total - sum(parts)))
    tel = {}
    for qn in ("t_RL", "r_L"):
        P = mpc(solvers["original_production"][qn])
        terms = [ref["R_dev"][qn] - ref["L"][qn], ref["R_code"][qn] - ref["R_dev"][qn], P - ref["R_code"][qn]]
        tel[f"original_{qn}"] = dict(terms_abs=[float(abs(x)) for x in terms], residual=tele(terms, P - ref["L"][qn]))
        Wd = mpc(solvers["WZ_delta_lead"][qn])
        terms = [ref["R_W"][qn] - ref["L"][qn], Wd - ref["R_W"][qn]]
        tel[f"WZ_delta_{qn}"] = dict(terms_abs=[float(abs(x)) for x in terms], residual=tele(terms, Wd - ref["L"][qn]))
        Wc = mpc(solvers["WZ_code_lead"][qn])
        terms = [ref["R_W"][qn] - ref["L"][qn], ref["R_W_code_lead"][qn] - ref["R_W"][qn], Wc - ref["R_W_code_lead"][qn]]
        tel[f"WZ_code_{qn}"] = dict(terms_abs=[float(abs(x)) for x in terms], residual=tele(terms, Wc - ref["L"][qn]))
    # signed phase decomposition (first order; phases add for small errors): term phases
    phase_terms = dict(
        original=dict(device=metrics(ref["R_dev"]["t_RL"], ref["L"]["t_RL"])["eps_phi"], lead=metrics(ref["R_code"]["t_RL"], ref["R_dev"]["t_RL"])["eps_phi"],
                      rounding=metrics(mpc(prod["t_RL"]), ref["R_code"]["t_RL"])["eps_phi"], total=metrics(mpc(prod["t_RL"]), ref["L"]["t_RL"])["eps_phi"]),
        WZ_code_lead=dict(device=metrics(ref["R_W"]["t_RL"], ref["L"]["t_RL"])["eps_phi"], lead=metrics(ref["R_W_code_lead"]["t_RL"], ref["R_W"]["t_RL"])["eps_phi"],
                          rounding=metrics(mpc(out_wz["WZ_code_lead"]["t_RL"]), ref["R_W_code_lead"]["t_RL"])["eps_phi"],
                          total=metrics(mpc(out_wz["WZ_code_lead"]["t_RL"]), ref["L"]["t_RL"])["eps_phi"]),
        WZ_delta_lead=dict(device=metrics(ref["R_W"]["t_RL"], ref["L"]["t_RL"])["eps_phi"],
                           rounding_and_lead_data=metrics(mpc(out_wz["WZ_delta_lead"]["t_RL"]), ref["R_W"]["t_RL"])["eps_phi"],
                           total=metrics(mpc(out_wz["WZ_delta_lead"]["t_RL"]), ref["L"]["t_RL"])["eps_phi"]))
    # ---- fixed-point trajectories (in W units; slices 1 .. N-1) ---------------------------------------------
    w_star_dev = 1 - 1 / lam_dev                                 # fixed point of the ORIGINAL represented device recurrence
    w_star_W = 1 - 1 / lam_W                                     # fixed point of the W/Z represented device recurrence
    traj = {}
    tr_rows = []
    w_orig = [1 + mpc(d) / mt for d in rep["d_trace"][:N - 1]]
    for tag, ws, wstar in (("original", w_orig, w_star_dev), ("WZ_delta_lead", [mpc(w) for w in out_wz["WZ_delta_lead"]["w"]], w_star_W),
                           ("WZ_code_lead", [mpc(w) for w in out_wz["WZ_code_lead"]["w"]], w_star_W)):
        e = [w - wstar for w in ws]
        ea = [float(abs(x)) for x in e]
        traj[tag] = dict(n=len(e), w_star_re=float(mp.re(wstar)), w_star_im=float(mp.im(wstar)), abs_e_first=ea[0], abs_e_last=ea[-1],
                         max_abs_e=max(ea), mean_e_re=float(mp.re(sum(e) / len(e))), mean_e_im=float(mp.im(sum(e) / len(e))),
                         rel_max_abs_e=max(ea) / float(abs(wstar)))
        for j in list(range(0, len(e), 50)) + [len(e) - 1]:
            tr_rows.append(dict(formulation=tag, slice=j + 1, e_re=float(mp.re(e[j])), e_im=float(mp.im(e[j])), abs_e=ea[j]))
    # ---- reproduction gate -------------------------------------------------------------------------------
    committed = {"L->R": 1.5424928001550615e-10, "R->L": 1.542470040583065e-10}
    gate = dict(prod_eps_phi_tRL=phase_terms["original"]["total"], committed_LR=committed["L->R"],
                diff_LR=abs(phase_terms["original"]["total"] - committed["L->R"]),
                prod_eps_phi_tLR=metrics(mpc(prod["t_LR"]), ref["L"]["t_LR"])["eps_phi"], committed_RL=committed["R->L"])
    gate["diff_RL"] = abs(gate["prod_eps_phi_tLR"] - committed["R->L"])
    gate["passed"] = bool(gate["diff_LR"] <= 1e-15 and gate["diff_RL"] <= 1e-15)
    replica_vs_prod = float(max(abs(complex(replica[q]) - complex(prod[q])) for q in ("t_RL", "t_LR", "r_L", "r_R")))
    lat_doc = {kk: dict(name=v["name"], device_diagonal=str(mp.nstr(v["a"], 30)), c_device=mp.nstr(v["c_dev"], 30),
                        lead_root=mp.nstr(v["lam"], 30), sigma=mp.nstr(v["sigma"], 30), gamma=mp.nstr(v["gamma"], 30),
                        c_lead_effective=mp.nstr(v["c_lead"], 30), c_mismatch=float(v["c_mismatch"]), gamma_minus_2tImlam=float(v["gamma_vs_2tImlam"]),
                        mutually_consistent=v["consistent"], note=v["note"]) for kk, v in lat.items()}
    res = dict(
        experiment="W/Z RGF worst-point diagnostic (single point; experimental; not a production change)",
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"), git_head=git("rev-parse", "HEAD"),
        git_dirty=bool(git("status", "--porcelain")), python=sys.version.split()[0], numpy=np.__version__, mpmath=mp.__version__, dps=DPS,
        predictions_file="PREDICTIONS.md", predictions_sha256="sha256:" + hashlib.sha256(open(os.path.join(HERE, "PREDICTIONS.md"), "rb").read()).hexdigest(),
        point=dict(config="F", N_y=63, dx=dx, L_d=3.0, N=N, E=E, E_repr=repr(E), channel=CHANNEL, E_perp=Ep, t=t,
                   delta_fl=delta_fl, delta_exact=float((mE - mEp) * mdx ** 2), delta_fl_is_exact=bool(mp.mpf(delta_fl) == (mE - mEp) * mdx ** 2),
                   sterbenz_condition=bool(Ep / 2 <= E <= 2 * Ep), theta=float(2 * mp.asin(mp.sqrt(((mE - mEp) * mdx ** 2) / 2))),
                   sin_theta=float(mp.sqrt(((mE - mEp) * mdx ** 2) * (2 - (mE - mEp) * mdx ** 2))), a_dev=a_dev,
                   production_offdiag_max=off_t, production_min_rcond=float(r["min_rcond"])),
        lattices=lat_doc, R_W_identical_to_L=bool(lat["R_W"]["a"] == lat["L"]["a"] and lat["R_W"]["sigma"] == lat["L"]["sigma"]),
        reference_closed_form_vs_recursion=xcheck, reference_vs_lambda_power=analytic_check,
        reproduction_gate=gate, replica_vs_production_max_abs=replica_vs_prod,
        comparisons=rows, reference_differences=refdiff, telescoping=tel, signed_phase_terms=phase_terms, trajectories=traj,
        W_max=dict(WZ_delta_lead=out_wz["WZ_delta_lead"]["W_max"], WZ_code_lead=out_wz["WZ_code_lead"]["W_max"]),
        Z_N=dict(WZ_delta_lead=[out_wz["WZ_delta_lead"]["Z"].real, out_wz["WZ_delta_lead"]["Z"].imag], WZ_code_lead=[out_wz["WZ_code_lead"]["Z"].real, out_wz["WZ_code_lead"]["Z"].imag],
                 abs_delta_lead=abs(out_wz["WZ_delta_lead"]["Z"]), abs_code_lead=abs(out_wz["WZ_code_lead"]["Z"])),
        b1_criterion_unchanged=1e-10, wall_time_s=None)
    res["wall_time_s"] = time.time() - t0
    with open(os.path.join(OUT, "worst_point.json"), "w") as f:
        json.dump(res, f, indent=2)
    keys = list(rows[0])
    with open(os.path.join(OUT, "worst_point_comparisons.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for rr in rows:
            w.writerow({k_: (repr(v) if isinstance(v, float) else v) for k_, v in rr.items()})
    with open(os.path.join(OUT, "worst_point_trajectory.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(tr_rows[0]))
        w.writeheader()
        for rr in tr_rows:
            w.writerow({k_: (repr(v) if isinstance(v, float) else v) for k_, v in rr.items()})
    print(json.dumps(dict(gate=gate, replica_vs_prod=replica_vs_prod, R_W_identical_to_L=res["R_W_identical_to_L"], xcheck=xcheck,
                          analytic=analytic_check, phase_terms=phase_terms, telescoping=tel, traj=traj, W_max=res["W_max"], Z=res["Z_N"],
                          wall=res["wall_time_s"]), indent=1, default=str))
    for rr in rows:
        print(rr["id"], rr["solver"], "vs", rr["reference"], f"eps_phi={rr['eps_phi_tRL']:.4e} eps_A={rr['eps_A_tRL']:.3e} |q|={rr['abs_q_tRL']:.3e} "
              f"|r|={rr['abs_r_solver']:.3e} U={rr['unitarity_defect_solver']:.3e} Q={rr['reciprocity_defect_solver']:.3e}")
    for rr in refdiff:
        print("REF", rr["lattice"], "-", rr["minus"], f"eps_phi={rr['eps_phi_tRL']:.4e} eps_A={rr['eps_A_tRL']:.3e} |r_lat|={rr['abs_r_lattice']:.3e} U_lat={rr['unitarity_defect_lattice']:.2e}")


if __name__ == "__main__":
    main()
