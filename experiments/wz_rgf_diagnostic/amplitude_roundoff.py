r"""
EXPERIMENTAL W/Z amplitude-roundoff investigation (diagnostic only; predictions: PREDICTIONS_amplitude_roundoff.md, committed first).

Point: B1 config F (N_y = 63, dx = 1/1024, N = 3483), E = 125.69231487269683, clean guide (V_j = 0 exactly).
Parts:
  B  full 63-channel W/Z through the matrix code path (Delta_j = (E I - diag(E_perp) - V_j) dx^2, all lead roots);
  C  per-stage amplitude/phase decomposition of t_RL for every open channel by an exact complex-log telescoping,
     for several linear-solve methods, on the 63x63 path and the 1x1 channel-5 path;
  D  N = 1 and N = 2 edge cases, all channels; equivalence of four high-precision formulations.
Reference for every W/Z run: its own represented lattice R_W,n (delta_n = Delta_fl[n, n] exactly, consistent exact leads),
per channel (the clean device is exactly diagonal). Production runs are compared with their own lattice R_code,n
(stored a_dev,n, production double sigma/gamma), and with L_n (exact E - E_perp) for context.
Reads stationary / validate_rgf_b1 / smatrix_spec read-only; writes only under results_amplitude_roundoff/.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import time

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from stationary import build_device, scattering_matrix  # noqa: E402  (read-only use)
from validate_rgf_b1 import device, spec_energies  # noqa: E402
import wz_rgf as wz  # noqa: E402
import wz_variants as wv  # noqa: E402

DPS = 110
E_TARGET = "125.69231487269683"
K5 = 4                                    # channel n = 5 (0-based index)
OUT = os.path.join(HERE, "results_amplitude_roundoff")
METHODS_63 = ("solve", "inv", "lu", "diagdiv", "refine")
METHODS_1 = ("solve", "inv", "lu", "diagdiv", "refine", "refine_exact")


def mpc(z):
    z = complex(z)
    return mp.mpc(z.real, z.imag)


def f(x):
    return float(x)


# ---------------------------------------------------------------- exact scalar lattices (dimensionless, Ghat = t G)
def exact_root(delta):
    c = 1 - delta
    return c + 1j * mp.sqrt(1 - c * c) if delta > 0 else c - mp.sqrt(c * c - 1)


def exact_hat(x, mu, N):
    """Ghat of M/t = tridiag(1, 2x, 1) + mu (e1 e1^T + eN eN^T) by the Chebyshev closed form."""
    x = mp.mpc(x)
    rho = x + mp.sqrt(x * x - 1)
    U = lambda k: (rho ** (k + 1) - rho ** (-(k + 1))) / (rho - 1 / rho)
    det = U(N) + 2 * mu * U(N - 1) + mu * mu * U(N - 2)
    return dict(G11=(U(N - 1) + mu * U(N - 2)) / det, GN1=(-1) ** (N + 1) / det)


def mp_wz(delta, lam, N):
    lm = lam - 1
    if N == 1:
        Zi = 1 / (2 * delta + 2 * lm)
        return dict(G11=Zi, GN1=Zi, G1N=Zi, GNN=Zi)
    W = 2 * delta + lm
    G11 = 1 / (W - 1)
    Gj1 = G1j = G11
    for _ in range(1, N - 1):
        W = 2 * delta + W / (1 - W)
        P = 1 / (1 - W)
        G11 = G11 - G1j * P * Gj1
        Gj1, G1j = P * Gj1, G1j * P
    Zi = 1 / (2 * delta + lm + W / (1 - W))
    return dict(G11=G11 + G1j * Zi * Gj1, GN1=-Zi * Gj1, G1N=-G1j * Zi, GNN=Zi)


def mp_orig(x, mu, N):
    """Original recursion of stationary.rgf_endpoint_blocks in units t = 1 (D = 2x, + mu at the ends, D -= g)."""
    g = GN1 = G1N = G11 = None
    for j in range(N):
        D = mp.mpc(2 * x)
        if j == 0:
            D += mu
        if j == N - 1:
            D += mu
        if j > 0:
            D -= g
        gn = 1 / D
        if j == 0:
            g = GN1 = G1N = G11 = gn
            continue
        G11 = G11 + G1N * gn * GN1
        G1N, GN1, g = -G1N * gn, -gn * GN1, gn
    return dict(G11=G11, GN1=GN1, G1N=G1N, GNN=g)


def mp_direct(x, mu, N):
    M = mp.matrix(N, N)
    for i in range(N):
        M[i, i] = 2 * x
        if i + 1 < N:
            M[i, i + 1] = M[i + 1, i] = 1
    M[0, 0] += mu
    M[N - 1, N - 1] += mu
    Gi = M ** -1
    return dict(G11=Gi[0, 0], GN1=Gi[N - 1, 0], G1N=Gi[0, N - 1], GNN=Gi[N - 1, N - 1])


def maxdiff(a, b):
    return max(abs(mp.mpc(a[k]) - mp.mpc(b[k])) for k in ("G11", "GN1", "G1N", "GNN"))


def full_ref(G):
    return dict(G11=G["G11"], GN1=G["GN1"], G1N=G["GN1"], GNN=G["G11"])


def amp_phase(num, ref):
    q = mpc(num) / ref
    return f(abs(q) - 1), f(mp.arg(q))


def UQ(S):
    I = np.eye(S.shape[0])
    U = max(np.max(np.abs(S.conj().T @ S - I)), np.max(np.abs(S @ S.conj().T - I)))
    return float(U), float(np.max(np.abs(S - S.T)))


def git(*a):
    try:
        return subprocess.run(["git", *a], capture_output=True, text=True, cwd=ROOT, timeout=20).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


# ---------------------------------------------------------------- Part C: exact telescoping of ln(t_d / t_ref)
def telescope(rec, k, lead, delta_mp, lam, s_exact, t_ref, t_d, N):
    """Complex-log decomposition for channel k (diagonal case). Real parts: amplitude; imaginary parts: phase."""
    f1 = mpc(rec["G11_init"][k])
    Wd = [mpc(w) for w in rec["W_diag"][:, k]]
    Pd = [mpc(p) for p in rec["P_diag"][:, k]]
    form = [mp.log(f1 * (Wd[0] - 1))] + [mp.log(Pd[j] * (1 - Wd[j + 1])) for j in range(len(Pd))]
    winp = [mp.log((1 / (Wd[0] - 1)) / (-lam))] + [mp.log((1 / (1 - Wd[j + 1])) / lam) for j in range(len(Pd))]
    prod = f1
    for p in Pd:
        prod *= p
    Gm = mpc(rec["G_Nm1_1"][k])
    Zd, Zid = mpc(rec["Z"][k, k]), mpc(rec["Zi"][k, k])
    GN1 = mpc(rec["GhatN1"][k, k])
    Wl = mpc(rec["W_last"][k])
    Zx = 2 * delta_mp + mpc(lead["lam_minus_1"][k]) + Wl / (1 - Wl)
    Zstar = lam - 1 / lam
    sq = mp.mpf(float(np.sqrt(lead["sin_theta"][k])))
    T = dict(
        P_formation=mp.fsum(form),
        W_recurrence=mp.fsum(winp),
        propagation_multiplication=mp.log(Gm / prod),
        Z_inputs=mp.log(Zstar / Zx),
        Z_formation=mp.log(Zx / Zd),
        Z_inversion=mp.log(Zid * Zd),
        GN1_multiplication=mp.log(GN1 / (-Zid * Gm)),
        lead_sin_theta=mp.log(sq * sq / s_exact),
        fisher_lee_multiplication=mp.log(mpc(t_d) / (2j * sq * sq * GN1)),
    )
    total = mp.log(mpc(t_d) / t_ref)
    resid = abs(mp.fsum(T.values()) - total)
    step = dict(
        form_amp=[f(mp.re(v)) for v in form], form_phase=[f(mp.im(v)) for v in form],
        w_amp=[f(mp.re(v)) for v in winp], w_phase=[f(mp.im(v)) for v in winp])
    total_amp = f(mp.re(total))
    out = dict(total_amp=total_amp, total_phase=f(mp.im(total)), residual=f(resid),
               terms_amp={a: f(mp.re(b)) for a, b in T.items()}, terms_phase={a: f(mp.im(b)) for a, b in T.items()},
               per_step={kk: dict(mean=float(np.mean(v)), min=float(np.min(v)), max=float(np.max(v)),
                                  n_distinct=int(len(set(v)))) for kk, v in step.items()},
               n_distinct_P=int(len(set(complex(p) for p in rec["P_diag"][:, k]))),
               n_distinct_W=int(len(set(complex(w) for w in rec["W_diag"][:, k]))),
               predicted_from_P_formation=f(mp.re(T["P_formation"])))
    out["P_formation_fraction_of_amp"] = (out["terms_amp"]["P_formation"] / total_amp) if total_amp != 0 else None
    out["W_recurrence_fraction_of_amp"] = (out["terms_amp"]["W_recurrence"] / total_amp) if total_amp != 0 else None
    return out


def run_wz(Delta, lead, method, refs, cond=False):
    rec = wv.wz_blocks_variant(Delta, lead["lam_minus_1"], method=method, record=True, cond=cond)
    S = wz.wz_fisher_lee(rec, lead["sin_theta"], lead["open"])
    return rec, S


def main():
    mp.mp.dps = DPS
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    dev = device(63, 1 / 1024, 3.0)
    Es = [e for e in spec_energies(dev.basis.E_perp) if repr(float(e)) == E_TARGET]
    assert len(Es) == 1
    E = float(Es[0])
    Ep = np.asarray(dev.basis.E_perp, dtype=float)
    dx, N, t, ny = float(dev.dx), int(dev.n_slices), float(dev.t_x), dev.basis.n_modes
    assert (N, ny, dx, t) == (3483, 63, 1 / 1024, 2.0 ** 19)
    V = np.asarray(dev.V)
    assert float(np.max(np.abs(V))) == 0.0
    # ---- Part B inputs: full Delta_j and all lead roots ---------------------------------------------------
    Delta = (E * np.eye(ny) - np.diag(Ep) - V) * (dx * dx)
    assert Delta.shape == (N, ny, ny) and np.all(Delta == Delta[0]) and float(np.max(np.abs(Delta[0] - np.diag(np.diag(Delta[0]))))) == 0.0
    dfl = np.diag(Delta[0]).copy()
    lead = wz.lead_from_delta(dfl)
    is_open = lead["open"]
    o_idx = list(np.where(is_open)[0])
    closed_idx = [n for n in range(ny) if not is_open[n]]
    assert o_idx == [0, 1, 2, 3, 4] and K5 in o_idx
    dmp = [mp.mpf(float(d)) for d in dfl]
    lam = [exact_root(d) for d in dmp]
    s_ex = [mp.sqrt(d * (2 - d)) if d > 0 else mp.mpf(0) for d in dmp]

    # ---- exact references for the long device (closed form) + recursion cross-checks ----------------------
    ref = {}
    xcheck = dict(closed_vs_mp_wz=0.0, closed_vs_mp_orig=0.0, wz_vs_orig=0.0)
    for n in range(ny):
        R = full_ref(exact_hat(dmp[n] - 1, lam[n], N))
        A = mp_wz(dmp[n], lam[n], N)
        B = mp_orig(dmp[n] - 1, lam[n], N)
        scale = max(abs(R["G11"]), abs(R["GN1"]))
        xcheck["closed_vs_mp_wz"] = max(xcheck["closed_vs_mp_wz"], f(maxdiff(R, A) / scale))
        xcheck["closed_vs_mp_orig"] = max(xcheck["closed_vs_mp_orig"], f(maxdiff(R, B) / scale))
        xcheck["wz_vs_orig"] = max(xcheck["wz_vs_orig"], f(maxdiff(A, B) / scale))
        ref[n] = R
    t_ref = {n: 2j * s_ex[n] * ref[n]["GN1"] for n in o_idx}
    r_ref = {n: 2j * s_ex[n] * ref[n]["G11"] - 1 for n in o_idx}
    lam_power = {n: f(abs(t_ref[n] - lam[n] ** (N - 1))) for n in o_idx}
    t_part0 = time.time() - t0

    # ---- closed-channel root handling --------------------------------------------------------------------------
    from smatrix_spec import lead_modes  # read-only use
    m = lead_modes(E, Ep, dx)
    c_code = 1.0 - (E - Ep) / (2.0 * t)                     # identical expression to lead_modes (production c)
    root = []
    for n in closed_idx:
        lm_ex = lam[n] - 1
        own = mp.mpf(float(c_code[n]))
        lam_own = own - mp.sqrt(own * own - 1)
        root.append(dict(channel=n + 1, delta=float(dfl[n]),
                         wz_rel_err_lam_minus_1=f(abs(mpc(lead["lam_minus_1"][n]) - lm_ex) / abs(lm_ex)),
                         prod_rel_err_lam_minus_1_vs_RW=f(abs(mpc(m["lam"][n]) - 1 - lm_ex) / abs(lm_ex)),
                         prod_rel_err_lam_minus_1_vs_own_c=f(abs(mpc(m["lam"][n]) - 1 - (lam_own - 1)) / abs(lam_own - 1)),
                         wz_sin_theta=float(lead["sin_theta"][n]), prod_gamma=float(m["gamma"][n]),
                         wz_lambda_real=bool(np.imag(lead["lam_minus_1"][n]) == 0.0), abs_lambda=f(abs(lam[n]))))
    root_summary = dict(
        n_closed=len(closed_idx),
        wz_max_rel_err=max(r["wz_rel_err_lam_minus_1"] for r in root),
        prod_max_rel_err_vs_own_c=max(r["prod_rel_err_lam_minus_1_vs_own_c"] for r in root),
        prod_max_rel_err_vs_RW=max(r["prod_rel_err_lam_minus_1_vs_RW"] for r in root),
        n6=next(r for r in root if r["channel"] == 6), n63=next(r for r in root if r["channel"] == 63),
        all_closed_sin_theta_zero=all(r["wz_sin_theta"] == 0.0 for r in root),
        all_closed_gamma_zero_prod=all(r["prod_gamma"] == 0.0 for r in root),
        all_closed_lambda_real_in_wz=all(r["wz_lambda_real"] for r in root))

    # ---- Part B/C: runs ------------------------------------------------------------------------------------------
    runs, rows, recs = {}, [], {}
    path_inputs = {"63x63": (Delta, lead, list(range(ny))),
                   "1x1_ch5": (np.ascontiguousarray(Delta[:, K5:K5 + 1, K5:K5 + 1]),
                               {k: v[K5:K5 + 1] for k, v in lead.items()}, [K5])}
    for path, methods in (("63x63", METHODS_63), ("1x1_ch5", METHODS_1)):
        D_, L_, chans = path_inputs[path]
        for meth in methods:
            ts = time.time()
            rec, S = run_wz(D_, L_, meth, ref, cond=(path == "63x63" and meth == "solve"))
            U, Q = UQ(S["S"])
            opn = [n for n in chans if is_open[n]]
            per = {}
            for kk, n in enumerate(opn):
                tRL, tLR, rL, rR = S["t_RL"][kk, kk], S["t_LR"][kk, kk], S["r_L"][kk, kk], S["r_R"][kk, kk]
                aR, pR = amp_phase(tRL, t_ref[n])
                aL, pL = amp_phase(tLR, t_ref[n])
                li = chans.index(n)
                tel = telescope(rec, li, L_, dmp[n], lam[n], s_ex[n], t_ref[n], tRL, N)
                res_exact = max(f(abs(1 - (1 - mpc(w)) * mpc(p))) for w, p in zip(rec["W_diag"][1:, li], rec["P_diag"][:, li]))
                per[n + 1] = dict(eps_A_tRL=aR, eps_phi_tRL=pR, eps_A_tLR=aL, eps_phi_tLR=pL,
                                  abs_r_L=float(abs(rL)), abs_r_R=float(abs(rR)), abs_dr_L=f(abs(mpc(rL) - r_ref[n])),
                                  t_RL=[float(tRL.real), float(tRL.imag)], res_P_exact_max=res_exact, telescoping=tel)
            # closed-channel Ghat blocks (63 path)
            closed_err = None
            if path == "63x63":
                ce = {}
                for n in closed_idx:
                    ce[n + 1] = {b: f(abs(mpc(rec["Ghat" + b[1:]][n, n]) - ref[n][b]) / abs(ref[n][b])) for b in ("G11", "GN1", "G1N", "GNN")}
                closed_err = dict(max_rel_G11=max(v["G11"] for v in ce.values()), max_rel_GNN=max(v["GNN"] for v in ce.values()),
                                  max_rel_GN1=max(v["GN1"] for v in ce.values()), max_rel_G1N=max(v["G1N"] for v in ce.values()),
                                  n6=ce[6], n63=ce[63])
            S_off = float(np.max(np.abs(S["S"] - np.diag(np.diag(S["S"])))))
            info = dict(path=path, method=meth, U_full_open_S=U, Q_full_open_S=Q, S_offdiag_max=S_off,
                        S_size=int(S["S"].shape[0]), offdiag_max_internal=rec["offdiag_max"],
                        res_P_double_max=float(np.max(rec["res_P"])), res_P_double_mean=float(np.mean(rec["res_P"])),
                        res_Z_double=rec["res_Z"], cond_Z=rec["cond_Z"], W_max=rec["W_max"],
                        abs_Z_ch5=float(abs(rec["Z"][chans.index(K5), chans.index(K5)])),
                        open_channels=per, closed_channels=closed_err, wall_s=time.time() - ts)
            if rec["cond_IW"].size:
                info["cond_IW_min"], info["cond_IW_max"] = float(rec["cond_IW"].min()), float(rec["cond_IW"].max())
            runs[f"{path}:{meth}"] = info
            recs[f"{path}:{meth}"] = rec
            print(f"{path:8s} {meth:12s} ch5 eps_A={per[5]['eps_A_tRL']:+.4e} eps_phi={per[5]['eps_phi_tRL']:+.4e} "
                  f"Pform={per[5]['telescoping']['terms_amp']['P_formation']:+.4e} U={U:.3e} Q={Q:.1e} ({time.time() - ts:.1f}s)", flush=True)

    # ---- faithfulness of the variant module and 1x1 vs 63x63 bitwise comparison ------------------------------------
    base63 = wz.wz_endpoint_blocks(Delta, lead["lam_minus_1"])
    base1 = wz.wz_endpoint_blocks(path_inputs["1x1_ch5"][0], path_inputs["1x1_ch5"][1]["lam_minus_1"])
    faithful = {
        "63x63": all(np.array_equal(base63[k], recs["63x63:solve"][k]) for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN")),
        "1x1_ch5": all(np.array_equal(base1[k], recs["1x1_ch5:solve"][k]) for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN"))}
    bitwise = {}
    for meth in METHODS_63:
        a, b = recs[f"63x63:{meth}"], recs[f"1x1_ch5:{meth}"]
        bitwise[meth] = dict(
            W_identical=bool(np.array_equal(a["W_diag"][:, K5], b["W_diag"][:, 0])),
            P_identical=bool(np.array_equal(a["P_diag"][:, K5], b["P_diag"][:, 0])),
            GhatN1_identical=bool(a["GhatN1"][K5, K5] == b["GhatN1"][0, 0]),
            Ghat11_identical=bool(a["Ghat11"][K5, K5] == b["Ghat11"][0, 0]),
            eps_A_63_minus_1=runs[f"63x63:{meth}"]["open_channels"][5]["eps_A_tRL"] - runs[f"1x1_ch5:{meth}"]["open_channels"][5]["eps_A_tRL"])
    distinct_P_ch5 = {}
    for key, rec in recs.items():
        li = 0 if key.startswith("1x1") else K5
        distinct_P_ch5[key] = [float(rec["P_diag"][0, li].real), float(rec["P_diag"][0, li].imag)]
    n_distinct_P_values = len({tuple(v) for v in distinct_P_ch5.values()})
    # correctly rounded P for channel 5 (from the double W of the solve run; W is constant) and its modulus defect
    w5 = mpc(recs["63x63:solve"]["W_diag"][1, K5])
    Pex = 1 / (1 - w5)
    P_cr = complex(float(mp.re(Pex)), float(mp.im(Pex)))
    P_cr_amp = f(abs(mpc(P_cr)) / abs(lam[K5]) - 1)
    P_cr_info = dict(P_correctly_rounded=[P_cr.real, P_cr.imag], amp_defect_correctly_rounded=P_cr_amp,
                     predicted_eps_A_if_all_P_correctly_rounded=(N - 1) * P_cr_amp,
                     which_methods_give_correctly_rounded_P={k: bool(complex(*v) == P_cr) for k, v in distinct_P_ch5.items()})

    # ---- production (original) on the same device: full S -----------------------------------------------------------
    rp = scattering_matrix(E, dev)
    mprod = rp["modes"]
    Up, Qp = UQ(rp["S"])
    prod_ch = {}
    for kk, n in enumerate(o_idx):
        a_dev = E - (2.0 * t + Ep[n])
        Rc = exact_hat(mp.mpf(a_dev) / (2 * t), -mpc(mprod["sigma"][n]) / t, N)
        tRc = 1j * mp.mpf(float(mprod["gamma"][n])) * Rc["GN1"] / t
        dL = (mp.mpf(E) - mp.mpf(float(Ep[n]))) * mp.mpf(dx) ** 2
        tL = exact_root(dL) ** (N - 1)
        aR, pR = amp_phase(rp["t_RL"][kk, kk], tRc)
        aL, pL = amp_phase(rp["t_RL"][kk, kk], tL)
        aW, pW = amp_phase(complex(runs["63x63:solve"]["open_channels"][n + 1]["t_RL"][0], runs["63x63:solve"]["open_channels"][n + 1]["t_RL"][1]), tL)
        prod_ch[n + 1] = dict(eps_A_vs_Rcode=aR, eps_phi_vs_Rcode=pR, eps_A_vs_L=aL, eps_phi_vs_L=pL,
                              wz_solve_eps_A_vs_L=aW, wz_solve_eps_phi_vs_L=pW,
                              RW_minus_L_phase=f(mp.arg(t_ref[n] / tL)), RW_minus_L_amp=f(abs(t_ref[n] / tL) - 1),
                              abs_r_L=float(abs(rp["r_L"][kk, kk])))
    production = dict(U_full_open_S=Up, Q_full_open_S=Qp, S_offdiag_max=float(np.max(np.abs(rp["S"] - np.diag(np.diag(rp["S"]))))),
                      min_rcond=float(rp["min_rcond"]), channels=prod_ch)

    # ---- Part D: N = 1 and N = 2 ----------------------------------------------------------------------------------------
    edge = {}
    for Ne in (1, 2):
        exact_equiv = 0.0
        refe = {}
        for n in range(ny):
            R = full_ref(exact_hat(dmp[n] - 1, lam[n], Ne))
            forms = (mp_wz(dmp[n], lam[n], Ne), mp_orig(dmp[n] - 1, lam[n], Ne), mp_direct(dmp[n] - 1, lam[n], Ne))
            sc = max(abs(R["G11"]), abs(R["GN1"]))
            exact_equiv = max(exact_equiv, max(f(maxdiff(R, F_) / sc) for F_ in forms), f(maxdiff(forms[0], forms[1]) / sc))
            refe[n] = R
        res = dict(exact_formulations_max_rel_diff=exact_equiv, methods={})
        for meth in METHODS_63:
            rec, S = run_wz(np.ascontiguousarray(Delta[:Ne]), lead, meth, refe)
            U, Q = UQ(S["S"])
            err = {}
            for n in range(ny):
                err[n] = max(f(abs(mpc(rec["Ghat" + b[1:]][n, n]) - refe[n][b]) / abs(refe[n][b])) for b in ("G11", "GN1", "G1N", "GNN"))
            topen = {}
            for kk, n in enumerate(o_idx):
                te = 2j * s_ex[n] * refe[n]["GN1"]
                a_, p_ = amp_phase(S["t_RL"][kk, kk], te)
                topen[n + 1] = dict(eps_A=a_, eps_phi=p_)
            res["methods"][meth] = dict(U_full_open_S=U, Q_full_open_S=Q, offdiag_max=rec["offdiag_max"],
                                        max_rel_err_open=max(err[n] for n in o_idx), rel_err_closed_n6=err[5], rel_err_closed_n63=err[62],
                                        max_rel_err_closed=max(err[n] for n in closed_idx), cond_Z=rec["cond_Z"], res_Z=rec["res_Z"],
                                        open_t=topen)
        if Ne == 2:
            pot0 = lambda X, Y: np.zeros_like(X)  # noqa: E731
            dev2 = build_device(pot0, 1.7, 1.7 + dx, dx, 1.0, ny)
            assert dev2.n_slices == 2
            r2 = scattering_matrix(E, dev2)
            U2, Q2 = UQ(r2["S"])
            pe = {}
            for kk, n in enumerate(o_idx):
                a_dev = E - (2.0 * t + Ep[n])
                Rc = exact_hat(mp.mpf(a_dev) / (2 * t), -mpc(r2["modes"]["sigma"][n]) / t, 2)
                gam = mp.mpf(float(r2["modes"]["gamma"][n]))
                te, re_ = 1j * gam * Rc["GN1"] / t, 1j * gam * Rc["G11"] / t - 1
                pe[n + 1] = dict(rel_err_t=f(abs(mpc(r2["t_RL"][kk, kk]) - te) / abs(te)),
                                 rel_err_r=f(abs(mpc(r2["r_L"][kk, kk]) - re_) / max(abs(re_), mp.mpf(10) ** -300)))
            res["production_N2"] = dict(U_full_open_S=U2, Q_full_open_S=Q2, channels=pe,
                                        max_rel_err_t=max(v["rel_err_t"] for v in pe.values()))
        edge[f"N={Ne}"] = res

    # ---- hypothesis H_A evaluation (pre-registered rule) -----------------------------------------------------------------
    verdicts = {}
    for key, info in runs.items():
        tl = info["open_channels"][5]["telescoping"]
        tot = tl["total_amp"]
        pf, wr = tl["terms_amp"]["P_formation"], tl["terms_amp"]["W_recurrence"]
        others = {k: v for k, v in tl["terms_amp"].items() if k != "P_formation"}
        supported = abs(tot - pf) <= 0.1 * abs(tot) + 2e-15 and abs(wr) < 0.1 * abs(pf)
        refuted = any(abs(v) >= 0.5 * abs(tot) for v in others.values()) and abs(tot) > 2e-15
        verdicts[key] = dict(total_amp=tot, P_formation=pf, W_recurrence=wr, supported=bool(supported), refuted=bool(refuted),
                             max_other_term=max(others.items(), key=lambda kv: abs(kv[1])))
    all_sup = all(v["supported"] for v in verdicts.values())
    any_ref = any(v["refuted"] for v in verdicts.values())
    H_A = "REFUTED" if any_ref else ("SUPPORTED" if all_sup else "PARTIAL")
    if n_distinct_P_values < 2 and H_A == "SUPPORTED":
        H_A = "SUPPORTED (consistent but not independently discriminated: all methods gave bit-identical P)"

    out = dict(
        experiment="W/Z amplitude-roundoff investigation (single point, clean; experimental; not a validation)",
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), git_head=git("rev-parse", "HEAD"),
        git_status_porcelain=git("status", "--porcelain"),
        python=sys.version.split()[0], numpy=np.__version__, mpmath=mp.__version__, dps=DPS,
        predictions_file="PREDICTIONS_amplitude_roundoff.md",
        predictions_sha256="sha256:" + hashlib.sha256(open(os.path.join(HERE, "PREDICTIONS_amplitude_roundoff.md"), "rb").read()).hexdigest(),
        point=dict(config="F", N_y=ny, dx=dx, N=N, E=E, t=t, n_open=len(o_idx), n_closed=len(closed_idx),
                   V_max_abs=0.0, Delta_offdiag_max=0.0, delta=[float(d) for d in dfl]),
        reference_crosschecks=dict(long_device=xcheck, t_ref_vs_lambda_power=lam_power),
        closed_root_handling=dict(summary=root_summary, per_channel=root),
        runs=runs, variant_faithful_to_wz_rgf=faithful, bitwise_1x1_vs_63x63=bitwise,
        ch5_first_P_by_run=distinct_P_ch5, n_distinct_P_values_ch5=n_distinct_P_values, correctly_rounded_P=P_cr_info,
        production_original=production, edge_cases=edge, H_A=dict(verdict=H_A, per_run=verdicts),
        b1_criterion_unchanged=1e-10, wall_time_s=time.time() - t0, reference_time_s=t_part0)
    with open(os.path.join(OUT, "amplitude_roundoff_results.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    # ---- CSV ------------------------------------------------------------------------------------------------------------
    cols = ["part", "N", "path", "method", "channel", "eps_A_tRL", "eps_phi_tRL", "eps_A_tLR", "eps_phi_tLR", "abs_r_L", "abs_dr_L",
            "U_full_open_S", "Q_full_open_S", "res_P_double_max", "res_P_exact_max", "step_amp_P_formation_mean", "step_phase_P_formation_mean",
            "step_amp_W_mean", "T_amp_P_formation", "T_amp_W_recurrence", "T_amp_propagation", "T_amp_Z_total", "T_amp_fisher_lee_total",
            "T_phase_P_formation", "T_phase_W_recurrence", "T_phase_propagation", "telescoping_residual"]
    with open(os.path.join(OUT, "amplitude_roundoff_comparisons.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for key, info in runs.items():
            for ch, v in info["open_channels"].items():
                tl = v["telescoping"]
                ta, tp = tl["terms_amp"], tl["terms_phase"]
                w.writerow(dict(part="B/C", N=N, path=info["path"], method=info["method"], channel=ch,
                                eps_A_tRL=v["eps_A_tRL"], eps_phi_tRL=v["eps_phi_tRL"], eps_A_tLR=v["eps_A_tLR"], eps_phi_tLR=v["eps_phi_tLR"],
                                abs_r_L=v["abs_r_L"], abs_dr_L=v["abs_dr_L"], U_full_open_S=info["U_full_open_S"], Q_full_open_S=info["Q_full_open_S"],
                                res_P_double_max=info["res_P_double_max"], res_P_exact_max=v["res_P_exact_max"],
                                step_amp_P_formation_mean=tl["per_step"]["form_amp"]["mean"], step_phase_P_formation_mean=tl["per_step"]["form_phase"]["mean"],
                                step_amp_W_mean=tl["per_step"]["w_amp"]["mean"],
                                T_amp_P_formation=ta["P_formation"], T_amp_W_recurrence=ta["W_recurrence"], T_amp_propagation=ta["propagation_multiplication"],
                                T_amp_Z_total=ta["Z_inputs"] + ta["Z_formation"] + ta["Z_inversion"] + ta["GN1_multiplication"],
                                T_amp_fisher_lee_total=ta["lead_sin_theta"] + ta["fisher_lee_multiplication"],
                                T_phase_P_formation=tp["P_formation"], T_phase_W_recurrence=tp["W_recurrence"],
                                T_phase_propagation=tp["propagation_multiplication"], telescoping_residual=tl["residual"]))
        for ch, v in production["channels"].items():
            w.writerow(dict(part="B (production original vs R_code)", N=N, path="63x63", method="production", channel=ch,
                            eps_A_tRL=v["eps_A_vs_Rcode"], eps_phi_tRL=v["eps_phi_vs_Rcode"], abs_r_L=v["abs_r_L"],
                            U_full_open_S=production["U_full_open_S"], Q_full_open_S=production["Q_full_open_S"]))
        for Ne, res in edge.items():
            for meth, v in res["methods"].items():
                for ch, tv in v["open_t"].items():
                    w.writerow(dict(part="D", N=Ne[2:], path="63x63", method=meth, channel=ch, eps_A_tRL=tv["eps_A"], eps_phi_tRL=tv["eps_phi"],
                                    U_full_open_S=v["U_full_open_S"], Q_full_open_S=v["Q_full_open_S"]))
    print("H_A:", H_A, "| faithful:", faithful, "| distinct P:", n_distinct_P_values, "| wall %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
