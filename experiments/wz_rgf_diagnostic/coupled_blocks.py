r"""
EXPERIMENTAL coupled-block (non-commuting) validation of the W/Z recursion: diagnostic only, NOT a production change.
Design and pre-registered rules: PREDICTIONS_coupled_blocks.md (design commit 7d692608f068f2cd80fbf2e6472d271e7c0134bc).

Solvers (double, unmodified, read-only):  A = wz_rgf.wz_endpoint_blocks (+ instrumented copy wz_mutants.wz_instrumented,
checked bit-identical);  B = stationary.scattering_matrix on a directly built stationary.Device.
References C: the full block-tridiagonal matrix M = (diagonal blocks A_j, off-diagonal +tI) minus Sigma_L on slice 1 and minus
Sigma_R on slice N (two separate operations), assembled in mpmath from the matrix DEFINITION for each represented lattice
(L, R_W, R_orig) and inverted directly; prefix systems = same routine without Sigma_R. No recursion code (W/Z, original or
mutant) is used to build any reference.

Usage:  python3 coupled_blocks.py                 (76 executions -> results_coupled_blocks/)
        python3 coupled_blocks.py --smoke OUTDIR  (smoke subset with explicit assertions, outputs to OUTDIR)
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
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

from stationary import Device, device_hamiltonian_blocks, scattering_matrix, transverse_basis  # noqa: E402  (read-only)
from smatrix_spec import lead_modes  # noqa: E402  (read-only)
import wz_rgf as wz  # noqa: E402  (read-only)
import wz_mutants as wzm  # noqa: E402

DESIGN_COMMIT = "7d692608f068f2cd80fbf2e6472d271e7c0134bc"
DPS, DPS_HI = 60, 100
EPS = 2.0 ** -52
C_SAFE = 10.0
FAIL_THR = 1e-10
dx, Ly, v0 = 1 / 8, 1.0, 4.0
t = 1.0 / (2.0 * dx * dx)
EP_DESIGN = [4.934802200544679, 19.739208802178716, 44.41321980490211]
C2 = dict(mu=(0.5, -0.25, 1, 0), alpha=(1, 2, -1.5, 0.5), beta=(3, -1, 2.5, 1.75), gamma=(1, 0.5, -2, 1.25))
C3 = dict(mu=(0.5, -0.25, 1, 0), a=(1, 2, -1.5, 0.5), a2=(0.25, -0.5, 0.75, 1), b12=(3, -1, 2.5, 1.75), b23=(-2, 1.5, 0.5, -1),
          b13=(1, 0, -1.25, 2), g12=(1, 0.5, -2, 1.25), g23=(-0.5, 1, 0.75, -1.5))
ENER = {"e2a": (2, 12.0), "e2b": (2, 30.0), "e3a": (3, 12.0), "e3b": (3, 30.0), "e3c": (3, 50.0)}
SCALES = (1.0, 2.0 ** -4, 2.0 ** -8, 2.0 ** -12)


# ----------------------------------------------------------------------------------------------------------- devices
def V(ny, j, cls, s=1.0):
    """Identical to the definition in PREDICTIONS_coupled_blocks.md section 2 / coupled_preflight.py."""
    k = j - 1
    if ny == 2:
        c = C2
        A = np.diag([1.0, -1.0]); B = np.array([[0, 1.0], [1.0, 0]]); C = np.array([[0, -1j], [1j, 0]])
        M = c["mu"][k] * np.eye(2) + c["alpha"][k] * A + 0j
        if cls in "RH":
            M = M + c["beta"][k] * B
        if cls == "H":
            M = M + c["gamma"][k] * C
    else:
        c = C3
        E = lambda a, b: np.eye(3)[:, [a]] @ np.eye(3)[[b], :]  # noqa: E731
        S = lambda a, b: E(a, b) + E(b, a)  # noqa: E731
        K = lambda a, b: 1j * (E(b, a) - E(a, b))  # noqa: E731
        M = c["mu"][k] * np.eye(3) + c["a"][k] * np.diag([1.0, 0, -1]) + c["a2"][k] * np.diag([1.0, -2, 1]) + 0j
        if cls in "RH":
            M = M + c["b12"][k] * S(0, 1) + c["b23"][k] * S(1, 2) + c["b13"][k] * S(0, 2)
        if cls == "H":
            M = M + c["g12"][k] * K(0, 1) + c["g23"][k] * K(1, 2)
    return v0 * s * M


def configs():
    ex = []
    for eid in ENER:
        for cls in "RHD":
            for N in (1, 2, 3, 4):
                ex.append(dict(matrix="main", eid=eid, cls=cls, N=N, s=1.0))
    for eid in ("e2b", "e3b"):
        for N in (3, 4):
            for s in SCALES:
                ex.append(dict(matrix="scaling", eid=eid, cls="R", N=N, s=s))
    return ex


def key(c):
    return f"{c['eid']}/{c['cls']}/N{c['N']}" + ("" if c["s"] == 1.0 else f"/s=2^{int(round(math.log2(c['s'])))}")


# ----------------------------------------------------------------------------------------------------------- mp helpers
def mpc(z):
    z = complex(z)
    return mp.mpc(z.real, z.imag)


def to_mp(A):
    A = np.asarray(A)
    return mp.matrix([[mpc(x) for x in row] for row in A])


def to_np(M):
    return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])


def fro(M):
    return mp.sqrt(mp.fsum(abs(M[i, j]) ** 2 for i in range(M.rows) for j in range(M.cols)))


def maxabs(M):
    return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))


def cond2(M):
    return float(np.linalg.cond(to_np(M)))


def norm2(M):
    return float(np.linalg.norm(to_np(M), 2))


def blk(G, a, b, ny):
    return G[a * ny:(a + 1) * ny, b * ny:(b + 1) * ny]


def comm_rel(A, B):
    d = fro(A) * fro(B)
    return mp.mpf(0) if d == 0 else fro(A * B - B * A) / d


def exact_root(d):
    c = 1 - d
    return c + 1j * mp.sqrt(1 - c * c) if d > 0 else c - mp.sqrt(c * c - 1)


# ----------------------------------------------------------------------------------------------------------- lattices
def build_lattice(name, inp):
    """Represented lattice in physical units (evaluated at the current mp precision). Returns A blocks, t, Sigma, Fisher-Lee data."""
    ny, N = inp["ny"], inp["N"]
    tt = mp.mpf(t)
    I = mp.eye(ny)
    if name == "L":
        mE, mdx = mp.mpf(inp["E"]), mp.mpf(dx)
        Ep = [mp.mpf(float(e)) for e in inp["Ep"]]
        A = [mE * I - (2 * tt * I + mp.diag(Ep) + to_mp(inp["V"][j])) for j in range(N)]
        lam = [exact_root((mE - e) * mdx ** 2) for e in Ep]
        sig = [-tt * l for l in lam]
        fl = dict(kind="s", s=[mp.im(l) if mp.im(l) > 0 else mp.mpf(0) for l in lam])
    elif name == "R_W":
        A = [tt * (2 * to_mp(inp["Delta"][j]) - 2 * I) for j in range(N)]
        lam = [exact_root(mp.mpf(float(d))) for d in inp["dlead"]]
        sig = [-tt * l for l in lam]
        fl = dict(kind="s", s=[mp.im(l) if mp.im(l) > 0 else mp.mpf(0) for l in lam])
    elif name == "R_orig":
        A = [to_mp(inp["A_orig"][j]) for j in range(N)]
        sig = [mpc(x) for x in inp["sigma_code"]]
        lam = [-x / tt for x in sig]
        fl = dict(kind="gamma", gamma=[mp.mpf(float(g)) for g in inp["gamma_code"]])
    else:
        raise ValueError(name)
    return dict(name=name, A=A, t=tt, sig=sig, lam=lam, fl=fl, ny=ny)


def assemble(A, tt, sigL, sigR, right_lead=True):
    """M = blocktridiag(A_j; +tI) - Sigma_L on slice 1 - Sigma_R on slice n (separate operations; for n = 1 both on one block)."""
    n, ny = len(A), A[0].rows
    M = mp.matrix(n * ny, n * ny)
    for j in range(n):
        for a in range(ny):
            for b in range(ny):
                M[j * ny + a, j * ny + b] = A[j][a, b]
        if j + 1 < n:
            for a in range(ny):
                M[j * ny + a, (j + 1) * ny + a] = tt
                M[(j + 1) * ny + a, j * ny + a] = tt
    for a in range(ny):
        M[a, a] -= sigL[a]
    if right_lead:
        for a in range(ny):
            M[(n - 1) * ny + a, (n - 1) * ny + a] -= sigR[a]
    return M


def S_from_G(G, lat, N):
    """Exact open-channel S from Ghat (t units) of a lattice; convention P."""
    ny = lat["ny"]
    fl = lat["fl"]
    if fl["kind"] == "s":
        o = [n for n in range(ny) if fl["s"][n] > 0]
        f = [mp.sqrt(fl["s"][n]) for n in o]
        pref = lambda a, b: 2j * f[a] * f[b]  # noqa: E731
    else:
        o = [n for n in range(ny) if fl["gamma"][n] > 0]
        f = [mp.sqrt(fl["gamma"][n]) for n in o]
        pref = lambda a, b: 1j * f[a] * f[b] / lat["t"]  # noqa: E731
    no = len(o)
    B = {k: blk(G, *ij, ny) for k, ij in (("11", (0, 0)), ("N1", (N - 1, 0)), ("1N", (0, N - 1)), ("NN", (N - 1, N - 1)))}
    S = mp.matrix(2 * no, 2 * no)
    for a in range(no):
        for b in range(no):
            S[a, b] = pref(a, b) * B["11"][o[a], o[b]] - (1 if a == b else 0)
            S[a, no + b] = pref(a, b) * B["1N"][o[a], o[b]]
            S[no + a, b] = pref(a, b) * B["N1"][o[a], o[b]]
            S[no + a, no + b] = pref(a, b) * B["NN"][o[a], o[b]] - (1 if a == b else 0)
    return S, o


def UQ_mp(S):
    n = S.rows
    I = mp.eye(n)
    U = max(maxabs(S.H * S - I), maxabs(S * S.H - I))
    return U, maxabs(S - S.T)


def UQ_np(S):
    I = np.eye(S.shape[0])
    return float(max(np.max(np.abs(S.conj().T @ S - I)), np.max(np.abs(S @ S.conj().T - I)))), float(np.max(np.abs(S - S.T)))


def reference(lat, N, prefixes=True, checks=True, inp=None):
    """Direct-inverse reference + reference checks (a)-(c), and prefix references."""
    tt, ny = lat["t"], lat["ny"]
    M = assemble(lat["A"], tt, lat["sig"], lat["sig"], True)
    Minv = M ** -1
    G = tt * Minv
    out = dict(M=M, G=G, cond_M=cond2(M), norm2_G=norm2(G))
    if checks:
        n = M.rows
        Mc = M.copy()
        cols = [mp.lu_solve(Mc, mp.matrix([[1 if i == k else 0] for i in range(n)])) for k in range(n)]
        Glu = mp.matrix(n, n)
        for k in range(n):
            for i in range(n):
                Glu[i, k] = cols[k][i]
        out["check_a_lu_vs_inverse"] = float(fro(tt * Glu - G) / fro(G))
        out["check_c_residual"] = float(maxabs(M * Minv - mp.eye(n)))
        with mp.workdps(DPS_HI):
            lat_hi = build_lattice(lat["name"], inp)
            Ghi = lat_hi["t"] * assemble(lat_hi["A"], lat_hi["t"], lat_hi["sig"], lat_hi["sig"], True) ** -1
            out["check_b_dps100"] = float(fro(Ghi - G) / fro(Ghi))
        if N == 1:   # (d): single block receives both lead terms exactly once
            expect = lat["A"][0] - mp.diag(lat["sig"]) - mp.diag(lat["sig"])
            out["check_d_N1_two_leads"] = bool(all(M[i, j] == expect[i, j] for i in range(ny) for j in range(ny)))
            once = lat["A"][0] - mp.diag(lat["sig"])
            out["check_d_N1_not_single_lead"] = bool(any(M[i, i] != once[i, i] for i in range(ny)))
    pre = {}
    if prefixes:
        I = mp.eye(ny)
        for j in range(1, N):
            Mj = assemble(lat["A"][:j], tt, lat["sig"], lat["sig"], right_lead=False)
            Gp = tt * Mj ** -1
            Gjj = blk(Gp, j - 1, j - 1, ny)
            d = dict(G=Gp, cond=cond2(Mj), norm2=norm2(Gp), W=I + Gjj ** -1, P=-Gjj,
                     G11=blk(Gp, 0, 0, ny), Gj1=blk(Gp, j - 1, 0, ny), G1j=blk(Gp, 0, j - 1, ny))
            if checks:
                with mp.workdps(DPS_HI):
                    lat_hi = build_lattice(lat["name"], inp)
                    Ghi = lat_hi["t"] * assemble(lat_hi["A"][:j], lat_hi["t"], lat_hi["sig"], lat_hi["sig"], right_lead=False) ** -1
                    d["check_b_dps100"] = float(fro(Ghi - Gp) / fro(Ghi))
                d["check_c_residual"] = float(maxabs(Mj * (Gp / tt) - mp.eye(Mj.rows)))
                Mc = Mj.copy()
                n = Mj.rows
                Glu = mp.matrix(n, n)
                for k in range(n):
                    col = mp.lu_solve(Mc, mp.matrix([[1 if i == k else 0] for i in range(n)]))
                    for i in range(n):
                        Glu[i, k] = col[i]
                d["check_a_lu_vs_inverse"] = float(fro(tt * Glu - Gp) / fro(Gp))
            pre[j] = d
    out["prefix"] = pre
    out["Z"] = blk(G, N - 1, N - 1, ny) ** -1
    return out


# ------------------------------------------------------------------------------------- supplementary mp recursions (E2)
def mp_wz(lat, N):
    """W/Z equations in mpmath on lattice R_W (Delta = A/(2t) + I; Lambda = -Sigma/t)."""
    ny, tt = lat["ny"], lat["t"]
    I = mp.eye(ny)
    Dl = [A / (2 * tt) + I for A in lat["A"]]
    Lm = mp.diag([-s / tt for s in lat["sig"]]) - I
    if N == 1:
        Zi = (2 * Dl[0] + 2 * Lm) ** -1
        return dict(G11=Zi, GN1=Zi, G1N=Zi, GNN=Zi)
    W = 2 * Dl[0] + Lm
    G11 = (W - I) ** -1
    Gj1 = G1j = G11
    for j in range(1, N - 1):
        W = 2 * Dl[j] + (I - W) ** -1 * W
        P = (I - W) ** -1
        G11 = G11 - G1j * P * Gj1
        Gj1, G1j = P * Gj1, G1j * P
    Zi = (2 * Dl[N - 1] + Lm + (I - W) ** -1 * W) ** -1
    return dict(G11=G11 + G1j * Zi * Gj1, GN1=-Zi * Gj1, G1N=-G1j * Zi, GNN=Zi)


def mp_orig(lat, N):
    """Original RGF recursion in mpmath on a lattice (D_j = A_j - Sigma at ends - t^2 g); returns Ghat = t G."""
    ny, tt = lat["ny"], lat["t"]
    Sg = mp.diag(lat["sig"])
    g = Gj1 = G1j = G11 = None
    for j in range(N):
        D = lat["A"][j].copy()
        if j == 0:
            D = D - Sg
        if j == N - 1:
            D = D - Sg
        if j > 0:
            D = D - tt * tt * g
        gn = D ** -1
        if j == 0:
            g = Gj1 = G1j = G11 = gn
            continue
        G11 = G11 + tt * tt * (G1j * gn * Gj1)
        G1j = -tt * (G1j * gn)
        Gj1 = -tt * (gn * Gj1)
        g = gn
    return {k: tt * v for k, v in dict(G11=G11, GN1=Gj1, G1N=G1j, GNN=g).items()}


# ----------------------------------------------------------------------------------------------------------- metrics
def endpoint_refs(G, N, ny):
    return dict(G11=blk(G, 0, 0, ny), GN1=blk(G, N - 1, 0, ny), G1N=blk(G, 0, N - 1, ny), GNN=blk(G, N - 1, N - 1, ny))


def block_err(X, Xr, scale_full):
    D = to_mp(X) - Xr
    nX = fro(Xr)
    mx = maxabs(Xr)
    floor = max(mx * mp.mpf("1e-6"), mp.mpf(10) ** (-300))
    ent = max(abs(D[i, j]) / max(abs(Xr[i, j]), floor) for i in range(D.rows) for j in range(D.cols))
    return dict(e_glob=float(fro(D) / scale_full), e_rel=float(fro(D) / nX) if nX > 0 else None, e_abs=float(maxabs(D)),
                e_entry=float(ent), ref_fro=float(nX))


def bounds(ref, N, ny, n_open):
    """Section 5 bounds from reference quantities."""
    I = mp.eye(ny)
    pre = ref["prefix"]
    condIW = {j: cond2(I - pre[j]["W"]) for j in pre}
    normP = {j: norm2(pre[j]["P"]) for j in pre}
    normW = {j: norm2(pre[j]["W"]) for j in pre}
    cZ = cond2(ref["Z"])
    K = max([cZ] + list(condIW.values()))                       # cond(W_1 - I) = cond(I - W_1) is included
    g = max([1.0] + list(normP.values())) ** 2
    n = N * ny
    tau_end = max(1e-13, C_SAFE * n * EPS * max(ref["cond_M"], K) * g)
    out = dict(K_Schur=K, g=g, cond_Z=cZ, cond_IW=condIW, normP=normP, tau_end=tau_end)
    tp, tw, tP = {}, {}, {}
    for j in pre:
        Kj = max(condIW[k] for k in pre if k <= j)
        gj = max([1.0] + [normP[k] for k in pre if k <= j]) ** 2
        tp[j] = max(1e-13, C_SAFE * j * ny * EPS * max(pre[j]["cond"], Kj) * gj)
        tw[j] = C_SAFE * ny * EPS * j * g ** (j - 1) * (1 + K) * max([1.0] + [normW[k] for k in pre if k <= j])
        tP[j] = normP[j] ** 2 * tw[j] + C_SAFE * ny * EPS * condIW[j] * normP[j]
    out.update(tau_pre=tp, tau_W=tw, tau_P=tP)
    return out


def tau_S(ref, b, smax):
    return max(1e-13, 2 * smax * ref["norm2_G"] * b["tau_end"]) + 4 * EPS


def rho(A, X):
    A = np.asarray(A); X = np.asarray(X)
    R = np.eye(A.shape[0]) - A @ X
    return float(np.linalg.norm(R) / (np.linalg.norm(A) * np.linalg.norm(X)))


def arrays_hash(arrs):
    h = hashlib.sha256()
    for a in arrs:
        a = np.ascontiguousarray(np.asarray(a, dtype=complex))
        h.update(a.tobytes())
    return h.hexdigest()


def trace_dev(name, X, ref, Wref, N, ny):
    """Deviation of a traced intermediate from its exact reference (W/P absolute Frobenius; blocks e_glob)."""
    pre = ref["prefix"]
    if name.startswith("W_"):
        j = int(name[2:])
        return float(fro(to_mp(X) - pre[j]["W"]))
    if name.startswith("P_"):
        j = int(name[2:])
        return float(fro(to_mp(X) - pre[j]["P"]))
    if name.startswith("G11^(") or name.startswith("G_"):
        j = int(name.split("^(")[1][:-1])
        which = "G11" if name.startswith("G11^") else ("Gj1" if name.startswith(f"G_{j}1") else "G1j")
        return float(fro(to_mp(X) - pre[j][which]) / pre[j]["norm2"])
    if name == "Z":
        return float(fro(to_mp(X) - ref["Z"]) / fro(ref["Z"]))
    if name == "Zi":
        return float(fro(to_mp(X) - blk(ref["G"], N - 1, N - 1, ny)) / ref["norm2_G"])
    er = endpoint_refs(ref["G"], N, ny)
    return float(fro(to_mp(X) - er[name]) / ref["norm2_G"])


# ----------------------------------------------------------------------------------------------------------- one execution
def run_one(c):
    mp.mp.dps = DPS
    ny, E = ENER[c["eid"]]
    N, cls, s = c["N"], c["cls"], c["s"]
    basis = transverse_basis(Ly, ny)
    Ep = np.asarray(basis.E_perp, dtype=float)
    assert list(Ep) == EP_DESIGN[:ny]
    Vs = np.array([V(ny, j, cls, s) for j in range(1, N + 1)])
    if cls in "RD":
        assert np.all(Vs.imag == 0)
        Vs = Vs.real.copy()
    assert np.allclose(Vs, np.conj(np.transpose(Vs, (0, 2, 1))), atol=0, rtol=0)        # exactly Hermitian
    # ---- A: W/Z
    Delta = (E * np.eye(ny) - np.diag(Ep) - Vs) * (dx * dx)
    dlead = (E - Ep) * (dx * dx)
    lead = wz.lead_from_delta(dlead)
    GA = wz.wz_endpoint_blocks(Delta, lead["lam_minus_1"])
    SA = wz.wz_fisher_lee(GA, lead["sin_theta"], lead["open"])
    inst = wzm.wz_instrumented(Delta, lead["lam_minus_1"])
    faithful = all(np.array_equal(GA[k], inst[k]) for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN"))
    # ---- B: production original RGF
    dev = Device(xs=dx * np.arange(N), dx=dx, basis=basis, V=Vs, V_lead_max=0.0)
    rB = scattering_matrix(E, dev, min_closed=0)
    H, tB = device_hamiltonian_blocks(dev)
    assert tB == t
    eye = np.eye(ny, dtype=complex)
    A_orig = np.array([E * eye - H[j] for j in range(N)])                 # production's own formation (rgf_endpoint_blocks)
    modes = lead_modes(E, Ep, dx)
    assert np.array_equal(modes["sigma"], rB["modes"]["sigma"]) and np.array_equal(modes["gamma"], rB["modes"]["gamma"])
    GB = {k: t * rB["G"][k] for k in ("G11", "GN1", "G1N", "GNN")}
    n_open = int(np.sum(lead["open"]))
    assert n_open == int(rB["n_open"])
    inp = dict(ny=ny, N=N, E=E, Ep=Ep, V=Vs, Delta=Delta, dlead=dlead, A_orig=A_orig,
               sigma_code=modes["sigma"], gamma_code=modes["gamma"])
    # ---- references (C), reference checks first
    lats = {nm: build_lattice(nm, inp) for nm in ("L", "R_W", "R_orig")}
    refs = {nm: reference(lats[nm], N, prefixes=(nm != "L"), checks=True, inp=inp) for nm in lats}
    rc = {}
    for nm, r in refs.items():
        d = dict(a=r["check_a_lu_vs_inverse"], b=r["check_b_dps100"], c=r["check_c_residual"])
        for j, p in r["prefix"].items():
            d[f"prefix{j}_a"], d[f"prefix{j}_b"], d[f"prefix{j}_c"] = p["check_a_lu_vs_inverse"], p["check_b_dps100"], p["check_c_residual"]
        if N == 1:
            d["d_two_leads"], d["d_not_single"] = r["check_d_N1_two_leads"], r["check_d_N1_not_single_lead"]
        rc[nm] = d
    ok_abc = all(v <= 1e-55 for d in rc.values() for k, v in d.items() if isinstance(v, float))
    ok_d = all(d.get("d_two_leads", True) and d.get("d_not_single", True) for d in rc.values())
    rep_W = max(float(fro(refs["R_W"]["G"] - refs["L"]["G"]) / refs["L"]["norm2_G"]), 0.0)
    rep_O = max(float(fro(refs["R_orig"]["G"] - refs["L"]["G"]) / refs["L"]["norm2_G"]), 0.0)
    rep_bound = 10 * refs["L"]["cond_M"] * EPS
    ok_e = rep_W <= rep_bound and rep_O <= rep_bound
    # ---- E2 / E3 (supplementary / algebraic)
    rW, rO = refs["R_W"], refs["R_orig"]
    erW, erO = endpoint_refs(rW["G"], N, ny), endpoint_refs(rO["G"], N, ny)
    A2 = mp_wz(lats["R_W"], N)
    B2 = mp_orig(lats["R_orig"], N)
    E2 = dict(mp_wz_vs_C_RW=float(max(fro(A2[k] - erW[k]) for k in erW) / rW["norm2_G"]),
              mp_orig_vs_C_Rorig=float(max(fro(B2[k] - erO[k]) for k in erO) / rO["norm2_G"]))
    I = mp.eye(ny)
    pre = rW["prefix"]
    e3 = mp.mpf(0)
    for j in range(2, N):
        e3 = max(e3, maxabs(pre[j]["P"] * pre[j - 1]["Gj1"] - pre[j]["Gj1"]), maxabs(pre[j - 1]["G1j"] * pre[j]["P"] - pre[j]["G1j"]),
                 maxabs(pre[j - 1]["G11"] - pre[j - 1]["G1j"] * pre[j]["P"] * pre[j - 1]["Gj1"] - pre[j]["G11"]))
        Wrec = 2 * (lats["R_W"]["A"][j - 1] / (2 * lats["R_W"]["t"]) + I) + (I - pre[j - 1]["W"]) ** -1 * pre[j - 1]["W"]
        e3 = max(e3, maxabs(Wrec - pre[j]["W"]))
    Lam = mp.diag([-x / lats["R_W"]["t"] for x in lats["R_W"]["sig"]])
    DlN = lats["R_W"]["A"][N - 1] / (2 * lats["R_W"]["t"]) + I
    if N >= 2:
        Zf = 2 * DlN + (Lam - I) + (I - pre[N - 1]["W"]) ** -1 * pre[N - 1]["W"]
        W1f = 2 * (lats["R_W"]["A"][0] / (2 * lats["R_W"]["t"]) + I) + Lam - I
        e3 = max(e3, maxabs(Zf - rW["Z"]), maxabs(W1f - pre[1]["W"]))
    else:
        e3 = max(e3, maxabs(2 * DlN + 2 * (Lam - I) - rW["Z"]))
    E3 = float(e3)
    # ---- non-commutativity (from references)
    Dl = [lats["R_W"]["A"][j] / (2 * lats["R_W"]["t"]) + I for j in range(N)]
    kap = dict(kappa_W=[float(comm_rel(pre[j]["W"], pre[j + 1]["W"])) for j in range(1, N - 1)],
               kappa_chain=[float(comm_rel(pre[j]["P"], pre[j - 1]["Gj1"])) for j in range(2, N)],
               kappa_Z=float(comm_rel(rW["Z"] ** -1, pre[N - 1]["Gj1"])) if N >= 2 else None,
               kappa_DeltaLambda_max=float(max(comm_rel(D, Lam) for D in Dl)))
    m5 = []
    for j in range(3, N):
        X = pre[j - 1]["G1j"] * pre[j]["P"] * pre[j - 1]["Gj1"]; Y = pre[j - 1]["Gj1"] * pre[j]["P"] * pre[j - 1]["G1j"]
        m5.append(float(fro(X - Y) / fro(X)))
    kap["cond_M5"] = m5
    if N >= 2:
        Zi = rW["Z"] ** -1
        X = pre[N - 1]["G1j"] * Zi * pre[N - 1]["Gj1"]; Y = pre[N - 1]["Gj1"] * Zi * pre[N - 1]["G1j"]
        kap["cond_M8"] = float(fro(X - Y) / fro(X))
    else:
        kap["cond_M8"] = None
    # ---- bounds
    bW = bounds(rW, N, ny, n_open)
    bO = bounds(rO, N, ny, n_open)
    # ---- A vs C(R_W): endpoints, S, U, Q, symmetry
    def solver_metrics(Gd, Sd, ref, lat, b, sym_applicable):
        er = endpoint_refs(ref["G"], N, ny)
        blocks = {k: block_err(Gd[k], er[k], ref["norm2_G"]) for k in er}
        for k, v in blocks.items():
            v["e_rel_bound"] = b["tau_end"] * ref["norm2_G"] / v["ref_fro"] if v["ref_fro"] > 0 else None
        Sr, o = S_from_G(ref["G"], lat, N)
        Ur, Qr = UQ_mp(Sr)
        dS = float(maxabs(to_mp(Sd) - Sr))
        U, Q = UQ_np(Sd)
        smax = float(max(lat["fl"]["s"][n] for n in o)) if lat["fl"]["kind"] == "s" else float(max(lat["fl"]["gamma"][n] for n in o) / (2 * t))
        tS = tau_S(ref, b, smax)
        sym = float(np.max(np.abs(Gd["G1N"] - Gd["GN1"].T)))
        res = dict(blocks=blocks, max_e_glob=max(v["e_glob"] for v in blocks.values()), dS_max=dS, U=U, Q=Q, U_ref=float(Ur),
                   Q_ref=float(Qr), tau_end=b["tau_end"], tau_S=tS, U_bound=float(Ur) + 2 * (2 * len(o)) * tS,
                   sym_residual=sym, sym_bound=2 * b["tau_end"] * ref["norm2_G"], n_open=len(o))
        res["Q_metric"] = abs(Q - float(Qr))
        res["Q_bound"] = (2 * tS) if cls == "H" else float(Qr) + 2 * tS
        res["Q_value_judged"] = res["Q_metric"] if cls == "H" else Q
        res["sym_applicable"] = sym_applicable
        return res
    GdA = dict(G11=GA["Ghat11"], GN1=GA["GhatN1"], G1N=GA["Ghat1N"], GNN=GA["GhatNN"])
    mA = solver_metrics(GdA, SA["S"], rW, lats["R_W"], bW, cls in "RD")
    mB = solver_metrics(GB, rB["S"], rO, lats["R_orig"], bO, cls in "RD")
    mB["min_rcond"] = float(rB["min_rcond"])
    # ---- intermediate audit (N6) and residuals
    tr = dict(inst["trace"])
    inter = []
    for nm, X in inst["trace"]:
        if nm.startswith("W_") or nm.startswith("P_"):
            j = int(nm[2:])
            v = trace_dev(nm, X, rW, None, N, ny)
            bnd = bW["tau_W"][j] if nm.startswith("W_") else bW["tau_P"][j]
            inter.append(dict(name=nm, kind="W" if nm[0] == "W" else "P", value=v, bound=bnd))
        elif nm.startswith("G11^(") or nm.startswith("G_"):
            j = int(nm.split("^(")[1][:-1])
            inter.append(dict(name=nm, kind="prefix", value=trace_dev(nm, X, rW, None, N, ny), bound=bW["tau_pre"][j]))
        elif nm == "Z":
            inter.append(dict(name=nm, kind="Z", value=trace_dev(nm, X, rW, None, N, ny), bound=None))
    resid = []
    if N == 1:
        resid.append(dict(name="Z", rho=rho(tr["Z"], tr["Zi"])))
    else:
        resid.append(dict(name="W_1 - I", rho=rho(tr["W_1"] - np.eye(ny), tr["G11^(1)"])))
        for j in range(2, N):
            resid.append(dict(name=f"I - W_{j}", rho=rho(np.eye(ny) - tr[f"W_{j}"], tr[f"P_{j}"])))
        resid.append(dict(name="Z", rho=rho(tr["Z"], tr["Zi"])))
    # ---- mutants
    muts = []
    corr_trace = inst["trace"]
    corr_end = [inst[k] for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN")]
    for m in wzm.MUTANTS:
        r = wzm.wz_instrumented(Delta, lead["lam_minus_1"], mutant=m)
        ends = dict(G11=r["Ghat11"], GN1=r["GhatN1"], G1N=r["Ghat1N"], GNN=r["GhatNN"])
        er = endpoint_refs(rW["G"], N, ny)
        Dm = max(float(fro(to_mp(ends[k]) - er[k]) / rW["norm2_G"]) for k in er)
        Dvc = max(float(np.linalg.norm(ends[k] - inst[{"G11": "Ghat11", "GN1": "GhatN1", "G1N": "Ghat1N", "GNN": "GhatNN"}[k]]))
                  for k in er) / rW["norm2_G"]
        bit_ident = all(np.array_equal(a, b) for a, b in zip([r[k] for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN")], corr_end))
        first_bit, first_dev = None, None
        cmap = dict(corr_trace)
        for nm, X in r["trace"]:
            if first_bit is None and (nm not in cmap or not np.array_equal(X, cmap[nm])):
                first_bit = nm
            if first_dev is None and trace_dev(nm, X, rW, None, N, ny) > 1e-8:
                first_dev = nm
        if not wzm.on_path(m, N):
            cat = "off_path"
        elif m == "M12":
            cat = "equivalent"
        elif (m == "M5" and N == 3) or (m == "M8" and N == 2):
            cat = "identical_by_algebra"
        elif m in wzm.ORDER_MUTANTS and cls == "D":
            cat = "order_invisible_D"
        else:
            cat = "detectable"
        Sm = wz.wz_fisher_lee(r, lead["sin_theta"], lead["open"])["S"]
        Um, Qm = UQ_np(Sm)
        Sr, _ = S_from_G(rW["G"], lats["R_W"], N)
        muts.append(dict(mutant=m, category=cat, D_mut=Dm, D_vs_correct=Dvc, bit_identical_to_correct=bit_ident,
                         first_bitwise_difference=first_bit, first_deviation_gt_1e8=first_dev, predicted_first=wzm.PREDICTED_FIRST[m],
                         dS_max=float(maxabs(to_mp(Sm) - Sr)), U=Um, Q=Qm,
                         sym_residual=float(np.max(np.abs(r["Ghat1N"] - r["GhatN1"].T))),
                         hash=arrays_hash([r[k] for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN")])))
    hashes = dict(A=arrays_hash([GA[k] for k in ("Ghat11", "GhatN1", "Ghat1N", "GhatNN")] + [SA["S"]]),
                  B=arrays_hash([rB["G"][k] for k in ("G11", "GN1", "G1N", "GNN")] + [rB["S"]]),
                  mutants=arrays_hash([np.frombuffer(bytes.fromhex(x["hash"]), dtype=np.uint8).astype(complex) for x in muts]))
    return dict(config=dict(c, key=key(c), ny=ny, E=E, n_open=n_open, n_closed=ny - n_open), faithful=faithful,
                reference_checks=rc, ref_ok_abc=ok_abc, ref_ok_d=ok_d, rep_diff=dict(R_W_minus_L=rep_W, R_orig_minus_L=rep_O, bound=rep_bound),
                ref_ok_e=ok_e, E2=E2, E3=E3, kappa=kap,
                cond=dict(cond_M=rW["cond_M"], cond_M_orig=rO["cond_M"], cond_prefix_max=max([p["cond"] for p in pre.values()], default=None),
                          cond_I_minus_W_max=max(bW["cond_IW"].values(), default=None), cond_Z=bW["cond_Z"],
                          normP_max=max(bW["normP"].values(), default=None), K_Schur=bW["K_Schur"], g=bW["g"]),
                A=mA, B=mB, intermediates=inter, residuals=resid, mutants=muts, hashes=hashes,
                shapes=dict(Delta=list(Delta.shape), M=[rW["M"].rows, rW["M"].cols], S_A=list(SA["S"].shape), S_B=list(rB["S"].shape)))


# ----------------------------------------------------------------------------------------------------------- evaluation
def evaluate(results, preflight):
    main = [r for r in results if r["config"]["matrix"] == "main"]
    scal = [r for r in results if r["config"]["matrix"] == "scaling"]
    issues = dict(invalid=[], fail=[], anomalous=[], design_insufficient=[], production=[], notes=[])
    # validity
    for r in results:
        k = r["config"]["key"] + f" [{r['config']['matrix']}]"
        if not r["ref_ok_abc"]:
            issues["invalid"].append(f"{k}: reference check (a)-(c) failed: {r['reference_checks']}")
        if not r["ref_ok_d"]:
            issues["invalid"].append(f"{k}: reference check (d) (N = 1 two-lead assembly) failed")
        if not r["ref_ok_e"]:
            issues["invalid"].append(f"{k}: reference check (e) representation difference above 10 cond eps: {r['rep_diff']}")
        if not r["faithful"]:
            issues["invalid"].append(f"{k}: instrumented copy not bit-identical to wz_endpoint_blocks")
    for r in main:
        c, kp = r["config"], r["kappa"]
        if c["cls"] in "RH":
            conds = []
            if c["N"] >= 3:
                conds += kp["kappa_chain"]
                conds.append(kp["cond_M8"])
            if c["N"] >= 2:
                conds.append(kp["kappa_Z"])
            if c["N"] >= 4:
                conds += kp["cond_M5"]
            if conds and min(conds) < 1e-3:
                issues["invalid"].append(f"{c['key']}: ordering condition below 1e-3 ({min(conds):.2e})")
    # W/Z FAIL / anomalies (A)
    for r in results:
        k = r["config"]["key"] + f" [{r['config']['matrix']}]"
        A = r["A"]
        judged = A["tau_end"] <= 1e-11
        if not judged:
            issues["notes"].append(f"{k}: tau_end {A['tau_end']:.2e} > 1e-11: insufficiently conditioned for the fixed threshold")
        for bn, v in A["blocks"].items():
            if judged and v["e_glob"] > FAIL_THR:
                issues["fail"].append(f"{k}: A endpoint {bn} e_glob {v['e_glob']:.2e}")
            elif v["e_glob"] > A["tau_end"]:
                issues["anomalous"].append(f"{k}: A endpoint {bn} e_glob {v['e_glob']:.2e} > tau_end {A['tau_end']:.2e}")
            if v["e_rel_bound"] is not None and v["e_rel"] > v["e_rel_bound"]:
                issues["anomalous"].append(f"{k}: A endpoint {bn} e_rel {v['e_rel']:.2e} > scaled bound {v['e_rel_bound']:.2e}")
        if judged and A["dS_max"] > FAIL_THR:
            issues["fail"].append(f"{k}: A dS_max {A['dS_max']:.2e}")
        elif A["dS_max"] > A["tau_S"]:
            issues["anomalous"].append(f"{k}: A dS_max {A['dS_max']:.2e} > tau_S {A['tau_S']:.2e}")
        if A["U"] > A["U_bound"]:
            issues["anomalous"].append(f"{k}: A U {A['U']:.2e} > bound {A['U_bound']:.2e}")
        if A["Q_value_judged"] > A["Q_bound"]:
            issues["anomalous"].append(f"{k}: A Q metric {A['Q_value_judged']:.2e} > bound {A['Q_bound']:.2e}")
        if A["sym_applicable"] and A["sym_residual"] > A["sym_bound"]:
            issues["anomalous"].append(f"{k}: A symmetry residual {A['sym_residual']:.2e} > {A['sym_bound']:.2e}")
        for it in r["intermediates"]:
            if it["bound"] is None:
                continue
            if judged and it["value"] > FAIL_THR:
                issues["fail"].append(f"{k}: A intermediate {it['name']} {it['value']:.2e}")
            elif it["value"] > it["bound"]:
                issues["anomalous"].append(f"{k}: A intermediate {it['name']} {it['value']:.2e} > bound {it['bound']:.2e}")
        for rs in r["residuals"]:
            if rs["rho"] > 100 * r["config"]["ny"] * EPS:
                issues["anomalous"].append(f"{k}: residual rho({rs['name']}) {rs['rho']:.2e} > 100 Ny eps")
        B = r["B"]
        for bn, v in B["blocks"].items():
            if v["e_glob"] > FAIL_THR:
                issues["production"].append(f"{k}: B endpoint {bn} e_glob {v['e_glob']:.2e} > 1e-10")
            elif v["e_glob"] > B["tau_end"]:
                issues["production"].append(f"{k}: B endpoint {bn} e_glob {v['e_glob']:.2e} > tau_end {B['tau_end']:.2e} (anomalous)")
        if B["dS_max"] > B["tau_S"]:
            issues["production"].append(f"{k}: B dS_max {B['dS_max']:.2e} > tau_S {B['tau_S']:.2e}")
        if B["U"] > B["U_bound"]:
            issues["production"].append(f"{k}: B U {B['U']:.2e} > bound {B['U_bound']:.2e}")
        if B["Q_value_judged"] > B["Q_bound"]:
            issues["production"].append(f"{k}: B Q metric {B['Q_value_judged']:.2e} > bound {B['Q_bound']:.2e}")
        if B["sym_applicable"] and B["sym_residual"] > B["sym_bound"]:
            issues["production"].append(f"{k}: B symmetry residual {B['sym_residual']:.2e} > {B['sym_bound']:.2e}")
    # mutants (N7, main matrix)
    mut_summary = {}
    for r in main:
        c = r["config"]
        for m in r["mutants"]:
            s = mut_summary.setdefault(m["mutant"], dict(detectable=0, detected=0, first_match=0, identical=0, identical_ok=0,
                                                         invisible_D=0, invisible_D_ok=0, equivalent=0, equivalent_ok=0, off_path=0,
                                                         off_path_ok=0, min_D_detectable=None, max_D_invisible=None, max_D_equiv=None,
                                                         Q_blind_R=0, Q_blind_R_cases=0))
            k = f"{c['key']}: {m['mutant']}"
            cat = m["category"]
            if cat == "detectable":
                s["detectable"] += 1
                det = m["D_mut"] >= 1e-6
                s["detected"] += det
                s["min_D_detectable"] = m["D_mut"] if s["min_D_detectable"] is None else min(s["min_D_detectable"], m["D_mut"])
                fm = m["first_deviation_gt_1e8"] == m["predicted_first"]
                s["first_match"] += fm
                if not det:
                    issues["design_insufficient"].append(f"{k}: predicted detectable, D_mut {m['D_mut']:.2e} < 1e-6")
                if not fm:
                    issues["notes"].append(f"{k}: first deviating quantity {m['first_deviation_gt_1e8']} != predicted {m['predicted_first']}")
                if c["cls"] == "R":
                    s["Q_blind_R_cases"] += 1
                    s["Q_blind_R"] += m["Q"] <= r["A"]["Q_bound"]
            elif cat == "identical_by_algebra":
                s["identical"] += 1
                s["identical_ok"] += m["bit_identical_to_correct"]
                if not m["bit_identical_to_correct"]:
                    issues["design_insufficient"].append(f"{k}: 'identical by algebra' but not bit-identical (D_vs_correct {m['D_vs_correct']:.2e})")
            elif cat == "order_invisible_D":
                s["invisible_D"] += 1
                ok = m["D_mut"] <= 1e-15
                s["invisible_D_ok"] += ok
                s["max_D_invisible"] = m["D_mut"] if s["max_D_invisible"] is None else max(s["max_D_invisible"], m["D_mut"])
                if m["D_mut"] > 1e-12:
                    issues["design_insufficient"].append(f"{k}: pure order mutant visible in class D (D_mut {m['D_mut']:.2e})")
                elif not ok:
                    issues["notes"].append(f"{k}: class-D order mutant D_mut {m['D_mut']:.2e} > 1e-15 (predicted <= 1e-15); D_vs_correct {m['D_vs_correct']:.2e}")
            elif cat == "equivalent":
                s["equivalent"] += 1
                ok = m["D_mut"] <= r["A"]["tau_end"]
                s["equivalent_ok"] += ok
                s["max_D_equiv"] = m["D_mut"] if s["max_D_equiv"] is None else max(s["max_D_equiv"], m["D_mut"])
                if not ok:
                    issues["design_insufficient"].append(f"{k}: M12 D_mut {m['D_mut']:.2e} > tau_end {r['A']['tau_end']:.2e}")
            else:
                s["off_path"] += 1
                s["off_path_ok"] += m["bit_identical_to_correct"]
                if not m["bit_identical_to_correct"]:
                    issues["invalid"].append(f"{k}: off-path mutant changed the output (implementation error)")
    # N8 scaling
    slopes = []
    groups = {}
    for r in scal:
        c = r["config"]
        groups.setdefault((c["eid"], c["N"]), []).append(r)
    for (eid, N), rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: -r["config"]["s"])
        ss = np.array([r["config"]["s"] for r in rs])
        corr = [r["A"]["max_e_glob"] for r in rs]
        taus = [r["A"]["tau_end"] for r in rs]
        corr_ok = all(e <= tb for e, tb in zip(corr, taus))
        mono = (all(np.diff(corr) < 0) or all(np.diff(corr) > 0)) and max(corr) / max(min(corr), 1e-300) > 10
        kchain = [min(r["kappa"]["kappa_chain"]) for r in rs]
        kslope = float(np.polyfit(np.log(ss), np.log(kchain), 1)[0])
        for m in wzm.ORDER_MUTANTS:
            vals = [next(x for x in r["mutants"] if x["mutant"] == m) for r in rs]
            if any(v["category"] != "detectable" for v in vals):
                continue
            D = np.array([v["D_mut"] for v in vals])
            sl = float(np.polyfit(np.log(ss), np.log(D), 1)[0])
            gg = all(d >= 100 * tb for d, tb in zip(D, taus))
            ok = 0.8 <= sl <= 2.2 and gg
            slopes.append(dict(base=f"{eid}/R/N{N}", mutant=m, D=list(map(float, D)), slope=sl, D_ge_100tau=bool(gg), ok=bool(ok)))
            if not ok:
                issues["design_insufficient"].append(f"{eid}/R/N{N} {m}: scaling slope {sl:.2f} or D>=100 tau ({gg}) outside prediction")
        slopes.append(dict(base=f"{eid}/R/N{N}", mutant="(correct code)", D=list(map(float, corr)), tau_end=taus, within_tau=bool(corr_ok),
                           monotone_trend=bool(mono), kappa_chain_slope=kslope))
        if not corr_ok or mono:
            issues["notes"].append(f"{eid}/R/N{N}: correct-code scaling check within_tau={corr_ok} monotone_trend={mono}")
    # N9 duplicates
    dups = []
    by_key = {}
    for r in main:
        by_key[r["config"]["key"]] = r
    for r in scal:
        if r["config"]["s"] == 1.0:
            m = by_key[r["config"]["key"]]
            same = r["hashes"] == m["hashes"]
            dups.append(dict(key=r["config"]["key"], bit_identical=bool(same)))
    N9 = all(d["bit_identical"] for d in dups) and len(dups) == 4
    # N4/N5 vs preflight
    pf = {x["case"]: x for x in preflight}
    n4 = []
    for r in results:
        if r["config"]["matrix"] == "scaling" and r["config"]["s"] == 1.0:
            continue
        p = pf[r["config"]["key"]]
        cmp = {}
        for kk, mine in (("cond_M", r["cond"]["cond_M"]), ("cond_prefix_max", r["cond"]["cond_prefix_max"]),
                         ("cond_I_minus_W_max", r["cond"]["cond_I_minus_W_max"]), ("cond_Z", r["cond"]["cond_Z"]),
                         ("normP_max", r["cond"]["normP_max"])):
            if p[kk] is None:
                cmp[kk] = 0.0 if mine is None else float("inf")
            else:
                cmp[kk] = abs(mine - p[kk]) / abs(p[kk])
        kd = 0.0
        for kk in ("kappa_W", "kappa_chain", "cond_M5"):
            for a, b in zip(r["kappa"][kk], p[kk]):
                kd = max(kd, 0.0 if a == b == 0 else abs(a - b) / max(abs(b), 1e-300))
        for kk in ("kappa_Z", "kappa_DeltaLambda_max"):
            a, b = r["kappa"][kk], p[kk]
            if a is not None:
                kd = max(kd, 0.0 if a == b == 0 else abs(a - b) / max(abs(b), 1e-300))
        a, b = r["kappa"]["cond_M8"], p.get("cond_M8")
        if a is not None:
            kd = max(kd, 0.0 if a == b == 0 else abs(a - b) / max(abs(b), 1e-300))
        n4.append(dict(key=r["config"]["key"], max_rel_cond_diff=max(cmp.values()), max_rel_kappa_diff=kd))
    N4 = all(x["max_rel_cond_diff"] <= 1e-6 and x["max_rel_kappa_diff"] <= 1e-6 for x in n4)
    RH = [r for r in main if r["config"]["cls"] in "RH"]
    D_ = [r for r in main if r["config"]["cls"] == "D"]
    mins = dict(kappa_chain=min(v for r in RH for v in r["kappa"]["kappa_chain"]),
                kappa_W=min(v for r in RH for v in r["kappa"]["kappa_W"]),
                kappa_Z=min(r["kappa"]["kappa_Z"] for r in RH if r["config"]["N"] >= 2),
                kappa_DL=min(r["kappa"]["kappa_DeltaLambda_max"] for r in RH),
                M5=min(v for r in RH for v in r["kappa"]["cond_M5"]),
                M8_N3plus=min(r["kappa"]["cond_M8"] for r in RH if r["config"]["N"] >= 3))
    M8_N2_zero = all(r["kappa"]["cond_M8"] == 0 for r in main if r["config"]["N"] == 2)
    D_zero = all(all(v == 0 for v in r["kappa"]["kappa_W"] + r["kappa"]["kappa_chain"] + r["kappa"]["cond_M5"]) and
                 (r["kappa"]["kappa_Z"] or 0) == 0 and r["kappa"]["kappa_DeltaLambda_max"] == 0 and (r["kappa"]["cond_M8"] or 0) == 0 for r in D_)
    kslopes = [x["kappa_chain_slope"] for x in slopes if x["mutant"] == "(correct code)"]
    N5 = (mins["kappa_chain"] >= 2.0e-2 and mins["kappa_W"] >= 3.6e-2 and mins["kappa_Z"] >= 7.0e-2 and mins["kappa_DL"] >= 9.9e-2
          and mins["M5"] >= 3.2e-2 and mins["M8_N3plus"] >= 2.8e-2 and M8_N2_zero and D_zero and all(abs(s - 1) <= 0.05 for s in kslopes))
    # N1 / N6 / predictions summaries
    N1 = not any("A endpoint" in x or "A dS_max" in x or "A U" in x or "A Q" in x or "A symmetry" in x
                 for x in issues["anomalous"] + issues["fail"])
    N6 = not any("A intermediate" in x for x in issues["anomalous"] + issues["fail"])
    E1 = all(r["ref_ok_abc"] and r["ref_ok_d"] for r in results)
    E2 = max(max(r["E2"].values()) for r in results)
    E3 = max(r["E3"] for r in results)
    if E2 > 1e-50:
        issues["notes"].append(f"E2 (supplementary) max {E2:.2e} > 1e-50")
    if E3 > 1e-50:
        issues["anomalous"].append(f"E3 identities max {E3:.2e} > 1e-50")
    if not N4:
        issues["notes"].append("N4: experiment-time conditioning/kappa differ from the preflight by > 1e-6 relative")
    if not N5:
        issues["notes"].append("N5: commutator minima / zeros / slopes differ from the preflight statement")
    if not N9:
        issues["fail"].append("N9: duplicated s = 1 executions not bit-identical (determinism)")
    invalid = bool(issues["invalid"])
    if invalid:
        decision = "INVALID"
    elif issues["fail"]:
        decision = "FAIL"
    elif issues["design_insufficient"]:
        decision = "DESIGN INSUFFICIENT"
    elif issues["anomalous"]:
        decision = "ANOMALOUS (not a pass)"
    elif N1 and N6 and N9:
        decision = "PASS (pre-registered scope only)"
    else:
        decision = "NOT PASSED (N1/N6/N9 not all satisfied)"
    return dict(decision=decision, issues=issues, mutant_summary=mut_summary, slopes=slopes, duplicates=dups,
                N1=N1, N4=N4, N5=N5, N6=N6, N9=N9, E1=E1, E2_max=E2, E3_max=E3, kappa_minima=mins, M8_N2_zero=M8_N2_zero,
                D_class_all_zero=D_zero, preflight_comparison=n4)


def git(*a):
    try:
        return subprocess.run(["git", *a], capture_output=True, text=True, cwd=ROOT, timeout=20).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def write_outputs(results, ev, outdir, meta):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "coupled_blocks_results.json"), "w") as fh:
        json.dump(dict(meta=meta, evaluation=ev, executions=results), fh, indent=1,
                  default=lambda o: o.item() if hasattr(o, "item") else str(o))
    cols = ["key", "matrix", "eid", "cls", "N", "s", "n_open", "row", "category", "max_e_glob", "e_glob_G11", "e_glob_GN1", "e_glob_G1N",
            "e_glob_GNN", "tau_end", "dS_max", "tau_S", "U", "U_bound", "Q", "Q_ref", "sym_residual", "D_mut", "D_vs_correct",
            "bit_identical_to_correct", "first_deviation", "predicted_first", "cond_M", "kappa_chain_min", "kappa_Z"]
    with open(os.path.join(outdir, "coupled_blocks_comparisons.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in results:
            c = r["config"]
            base = dict(key=c["key"], matrix=c["matrix"], eid=c["eid"], cls=c["cls"], N=c["N"], s=c["s"], n_open=c["n_open"],
                        cond_M=r["cond"]["cond_M"], kappa_chain_min=min(r["kappa"]["kappa_chain"], default=""), kappa_Z=r["kappa"]["kappa_Z"])
            for sv in ("A", "B"):
                X = r[sv]
                w.writerow(dict(base, row=f"solver {sv}", max_e_glob=X["max_e_glob"], **{f"e_glob_{k}": v["e_glob"] for k, v in X["blocks"].items()},
                                tau_end=X["tau_end"], dS_max=X["dS_max"], tau_S=X["tau_S"], U=X["U"], U_bound=X["U_bound"], Q=X["Q"],
                                Q_ref=X["Q_ref"], sym_residual=X["sym_residual"]))
            for m in r["mutants"]:
                w.writerow(dict(base, row=f"mutant {m['mutant']}", category=m["category"], D_mut=m["D_mut"], D_vs_correct=m["D_vs_correct"],
                                bit_identical_to_correct=m["bit_identical_to_correct"], first_deviation=m["first_deviation_gt_1e8"],
                                predicted_first=m["predicted_first"], dS_max=m["dS_max"], U=m["U"], Q=m["Q"], sym_residual=m["sym_residual"]))
    write_report(results, ev, outdir, meta)


def fmt(x, p="%.1e"):
    return "—" if x is None else (p % x)


def write_report(results, ev, outdir, meta):
    main = [r for r in results if r["config"]["matrix"] == "main"]
    L = []
    a = L.append
    a("# Coupled-block (non-commuting) W/Z diagnostic — results (experimental; generated by coupled_blocks.py)\n")
    a(f"Design: `PREDICTIONS_coupled_blocks.md` (design commit `{DESIGN_COMMIT}`, sha256 `{meta['predictions_sha256']}`). "
      f"Run at HEAD `{meta['git_head']}`; generated {meta['generated_utc']}; wall {meta['wall_s']:.0f} s.\n")
    a("**B1 configuration F still FAILS the unchanged 1e-10 criterion. W/Z is NOT validated by this experiment, whatever its outcome:** "
      "the pre-registered claim is limited to these small deterministic coupled configurations.\n")
    a(f"## Decision (pre-registered rules, section 8): **{ev['decision']}**\n")
    a("| flag | value |\n|---|---|")
    for k in ("E1", "N1", "N4", "N5", "N6", "N9", "M8_N2_zero", "D_class_all_zero"):
        a(f"| {k} | {ev[k]} |")
    a(f"| E2 (supplementary) max | {ev['E2_max']:.1e} |\n| E3 max | {ev['E3_max']:.1e} |")
    a(f"| invalid / fail / anomalous / design-insufficient / production / notes | {len(ev['issues']['invalid'])} / {len(ev['issues']['fail'])} / "
      f"{len(ev['issues']['anomalous'])} / {len(ev['issues']['design_insufficient'])} / {len(ev['issues']['production'])} / {len(ev['issues']['notes'])} |\n")
    a("## 1. Configuration counts\n")
    nu = len({r["config"]["key"] for r in results})
    a(f"Executions: **{len(results)}** (main {len(main)}, scaling {len(results) - len(main)}); unique configurations: **{nu}**. "
      f"Duplicated s = 1 executions: {len(ev['duplicates'])}, bit-identical: {sum(d['bit_identical'] for d in ev['duplicates'])}.\n")
    a("## 2. Reference validation (run before any solver comparison)\n")
    mx = lambda f: max(f(r) for r in results)  # noqa: E731
    a("| check | worst over all executions, lattices, prefixes | threshold |\n|---|---|---|")
    a(f"| (a) mp.lu_solve vs mp.inverse | {mx(lambda r: max(v for d in r['reference_checks'].values() for k, v in d.items() if k.endswith('a') and isinstance(v, float))):.1e} | 1e-55 |")
    a(f"| (b) dps 100 vs dps 60 | {mx(lambda r: max(v for d in r['reference_checks'].values() for k, v in d.items() if k.endswith('b') and isinstance(v, float))):.1e} | 1e-55 |")
    a(f"| (c) mp residual | {mx(lambda r: max(v for d in r['reference_checks'].values() for k, v in d.items() if k.endswith('c') and isinstance(v, float))):.1e} | 1e-55 |")
    n1 = [r for r in results if r["config"]["N"] == 1]
    a(f"| (d) N = 1 block = A_1 − Σ_L − Σ_R exactly, and ≠ A_1 − Σ (single lead) | {sum(r['ref_ok_d'] for r in n1)}/{len(n1)} executions × 3 lattices | exact |")
    a(f"| (e) C(R_W) − C(L) | {mx(lambda r: r['rep_diff']['R_W_minus_L']):.1e} | 10·cond₂(M)·ε (≤ {mx(lambda r: r['rep_diff']['bound']):.1e}) |")
    a(f"| (e) C(R_orig) − C(L) | {mx(lambda r: r['rep_diff']['R_orig_minus_L']):.1e} | same |")
    a(f"| instrumented copy bit-identical to `wz_endpoint_blocks` (4 endpoint blocks) | {sum(r['faithful'] for r in results)}/{len(results)} | all |")
    a(f"| E2: mp W/Z on R_W, mp original on R_orig vs C (supplementary) | {ev['E2_max']:.1e} | 1e-50 |")
    a(f"| E3: prefix, W-recurrence, closing and N = 1 identities on the references | {ev['E3_max']:.1e} | 1e-50 |")
    a(f"| N4: conditioning / κ vs committed preflight (max rel. diff.) | cond {max(x['max_rel_cond_diff'] for x in ev['preflight_comparison']):.1e}, κ {max(x['max_rel_kappa_diff'] for x in ev['preflight_comparison']):.1e} | 1e-6 |\n")
    a("## 3. Solver A (W/Z) vs C(R_W) and solver B (production original) vs C(R_orig)\n")
    a("e_glob = ‖ΔX‖_F/‖Ĝ_full‖₂ (max over the four endpoint blocks); ΔS = max |S − S_ref|; U, Q of the full open-channel S "
      "(class H: Q judged as |Q − Q_ref|); sym = ‖Ĝ_1N − Ĝ_N1ᵀ‖_max (classes R, D).\n")
    a("| execution | open | cond₂(M) | A e_glob | A τ_end | A ΔS | A U | A Q (judged) | A sym | B e_glob | B τ_end | B ΔS | B U | B Q (judged) | B sym |")
    a("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in results:
        A, B = r["A"], r["B"]
        tag = r["config"]["key"] + (" (dup)" if r["config"]["matrix"] == "scaling" and r["config"]["s"] == 1.0 else "")
        a(f"| {tag} | {r['config']['n_open']} | {r['cond']['cond_M']:.1f} | {A['max_e_glob']:.1e} | {A['tau_end']:.1e} | {A['dS_max']:.1e} | {A['U']:.1e} | "
          f"{A['Q_value_judged']:.1e} | {A['sym_residual']:.1e} | {B['max_e_glob']:.1e} | {B['tau_end']:.1e} | {B['dS_max']:.1e} | {B['U']:.1e} | "
          f"{B['Q_value_judged']:.1e} | {B['sym_residual']:.1e} |")
    a("")
    allA = [v["e_glob"] for r in results for v in r["A"]["blocks"].values()]
    allB = [v["e_glob"] for r in results for v in r["B"]["blocks"].values()]
    a(f"Summary: A endpoint e_glob max {max(allA):.1e} (fraction in the predicted 1e-16..1e-14 band: {np.mean([(1e-16 <= x <= 1e-14) for x in allA]):.2f}; "
      f"exact zeros: {sum(x == 0 for x in allA)}); B max {max(allB):.1e}. Max own-relative block error A "
      f"{max(v['e_rel'] for r in results for v in r['A']['blocks'].values()):.1e}.\n")
    a("## 4. Intermediate audit of A (prefix blocks, W_j, P_j) and solve residuals\n")
    kinds = {}
    for r in results:
        for it in r["intermediates"]:
            if it["bound"] is None:
                continue
            d = kinds.setdefault(it["kind"], dict(n=0, max=0.0, max_ratio=0.0))
            d["n"] += 1
            d["max"] = max(d["max"], it["value"])
            d["max_ratio"] = max(d["max_ratio"], it["value"] / it["bound"])
    a("| quantity | n | max error | max error / own bound |\n|---|---|---|---|")
    for k, d in kinds.items():
        a(f"| {k} | {d['n']} | {d['max']:.1e} | {d['max_ratio']:.2f} |")
    zs = [it["value"] for r in results for it in r["intermediates"] if it["kind"] == "Z"]
    a(f"| Z_N (relative, no bound) | {len(zs)} | {max(zs):.1e} | — |")
    rr = [x["rho"] for r in results for x in r["residuals"]]
    a(f"\nSolve residuals ρ: max {max(rr):.1e} (anomaly threshold 100·N_y·ε ≥ {100 * 2 * EPS:.1e}).\n")
    a("## 5. Mutants (main matrix, 60 executions)\n")
    a("Categories: detectable (must give D_mut ≥ 1e-6 and the predicted first deviating quantity); identical_by_algebra (must be "
      "bit-identical to the correct code; never counted as a detection); order_invisible_D (pure order mutant in class D, "
      "predicted D_mut ≤ 1e-15; >1e-12 = design insufficient); equivalent (M12, D_mut ≤ τ_end); off_path (code not executed; "
      "must be bit-identical). D_mut is measured against C(R_W); D_vs_correct against the correct double run.\n")
    a("| mutant | detectable: detected / n | first-deviation match | min D_mut (detectable) | Q-blind in class R (detectable cases) | identical-by-algebra: bit-identical / n | class D invisible: ≤1e-15 / n (max D_mut) | M12 ≤ τ_end / n (max) | off-path bit-identical / n |")
    a("|---|---|---|---|---|---|---|---|---|")
    for m, s in ev["mutant_summary"].items():
        a(f"| {m} | {s['detected']}/{s['detectable']} | {s['first_match']}/{s['detectable']} | {fmt(s['min_D_detectable'])} | "
          f"{s['Q_blind_R']}/{s['Q_blind_R_cases']} | {s['identical_ok']}/{s['identical']} | {s['invisible_D_ok']}/{s['invisible_D']} ({fmt(s['max_D_invisible'])}) | "
          f"{s['equivalent_ok']}/{s['equivalent']} ({fmt(s['max_D_equiv'])}) | {s['off_path_ok']}/{s['off_path']} |")
    a("\n'Q-blind' counts detectable cases in class R whose mutant Q stays within the correct code's Q bound, i.e. reciprocity alone "
      "would not have flagged them.\n")
    a("## 6. Coupling scaling (N8)\n")
    a("| base | mutant | D_mut at s = 1, 2⁻⁴, 2⁻⁸, 2⁻¹² | log-log slope | all D ≥ 100 τ_end | ok |\n|---|---|---|---|---|---|")
    for x in ev["slopes"]:
        if x["mutant"] == "(correct code)":
            a(f"| {x['base']} | correct code (max e_glob) | {', '.join('%.1e' % v for v in x['D'])} | κ_chain slope {x['kappa_chain_slope']:.3f} | "
              f"within τ: {x['within_tau']} | monotone trend: {x['monotone_trend']} |")
        else:
            a(f"| {x['base']} | {x['mutant']} | {', '.join('%.1e' % v for v in x['D'])} | {x['slope']:.2f} | {x['D_ge_100tau']} | {x['ok']} |")
    a("\nOperationalization fixed before the run: 'no monotone trend' = not (strictly monotone over the 4 scales with max/min > 10).\n")
    a("## 7. Non-commutativity (from the references)\n")
    km = ev["kappa_minima"]
    a(f"Main R/H minima: κ_chain {km['kappa_chain']:.2e}, κ_W {km['kappa_W']:.2e}, κ_Z {km['kappa_Z']:.2e}, κ_[Δ,Λ] {km['kappa_DL']:.2e}, "
      f"M5 condition {km['M5']:.2e}, M8 condition (N ≥ 3) {km['M8_N3plus']:.2e}; M8 condition exactly 0 at N = 2: {ev['M8_N2_zero']}; "
      f"class D all exactly 0: {ev['D_class_all_zero']}.\n")
    a("## 8. Every invalidity, failure, anomaly, design insufficiency, production finding and note\n")
    for k in ("invalid", "fail", "anomalous", "design_insufficient", "production", "notes"):
        a(f"### {k} ({len(ev['issues'][k])})\n")
        for x in ev["issues"][k]:
            a(f"- {x}")
        if not ev["issues"][k]:
            a("- none")
        a("")
    with open(os.path.join(outdir, "coupled_blocks_report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def main():
    smoke = "--smoke" in sys.argv
    outdir = sys.argv[sys.argv.index("--smoke") + 1] if smoke else os.path.join(HERE, "results_coupled_blocks")
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    exs = configs()
    assert len(exs) == 76 and len({key(c) for c in exs}) == 72
    if smoke:
        want = {"e2a/R/N1", "e2a/R/N2", "e2a/R/N3", "e2a/R/N4", "e3b/H/N1", "e3b/H/N4", "e3a/D/N3"}
        exs = [c for c in exs if c["matrix"] == "main" and key(c) in want]
    # preflight (committed script, unchanged) -> consistency predictions N4/N5
    pf_path = os.path.join(outdir, "preflight.json")
    subprocess.run([sys.executable, "-I", os.path.join(HERE, "coupled_preflight.py"), pf_path], check=True, capture_output=True)
    preflight = json.load(open(pf_path))
    results = []
    for c in exs:
        ts = time.time()
        r = run_one(c)
        r["wall_s"] = time.time() - ts
        results.append(r)
        print(f"{r['config']['key']:22s} {c['matrix']:7s} faithful={r['faithful']} refABC={r['ref_ok_abc']} refD={r['ref_ok_d']} refE={r['ref_ok_e']} "
              f"A={r['A']['max_e_glob']:.1e} B={r['B']['max_e_glob']:.1e} ({r['wall_s']:.1f}s)", flush=True)
        if not (r["ref_ok_abc"] and r["ref_ok_d"] and r["ref_ok_e"]):
            print("REFERENCE CHECK FAILED -- stopping (no repair).", flush=True)
            json.dump(r, open(os.path.join(outdir, "reference_failure_raw.json"), "w"), indent=1, default=str)
            sys.exit(2)
    if smoke:
        for r in results:
            ny, N = r["config"]["ny"], r["config"]["N"]
            assert r["shapes"]["Delta"] == [N, ny, ny] and r["shapes"]["M"] == [N * ny, N * ny]
            assert r["shapes"]["S_A"] == r["shapes"]["S_B"] == [2 * r["config"]["n_open"]] * 2
            assert r["faithful"]
            if N == 1:
                assert r["ref_ok_d"]
        print("SMOKE: shapes, N=1 two-lead assembly, faithfulness, reference checks OK")
        json.dump(results, open(os.path.join(outdir, "smoke_results.json"), "w"), indent=1, default=str)
        return
    ev = evaluate(results, preflight)
    meta = dict(design_commit=DESIGN_COMMIT, git_head=git("rev-parse", "HEAD"), git_status=git("status", "--porcelain"),
                predictions_sha256=hashlib.sha256(open(os.path.join(HERE, "PREDICTIONS_coupled_blocks.md"), "rb").read()).hexdigest(),
                preflight_script_sha256=hashlib.sha256(open(os.path.join(HERE, "coupled_preflight.py"), "rb").read()).hexdigest(),
                python=sys.version.split()[0], numpy=np.__version__, mpmath=mp.__version__, dps=DPS, dps_hi=DPS_HI,
                generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), wall_s=time.time() - t0,
                n_executions=len(results), n_unique=len({r["config"]["key"] for r in results}))
    write_outputs(results, ev, outdir, meta)
    print("DECISION:", ev["decision"])
    for k, v in ev["issues"].items():
        print(f"  {k}: {len(v)}")


if __name__ == "__main__":
    main()
