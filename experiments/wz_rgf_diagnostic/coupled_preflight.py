"""Design preflight for PREDICTIONS_coupled_blocks.md (NOT the experiment).
Uses only the defined matrices/energies; exact quantities from direct mpmath inverses of the full and prefix matrices of the
represented W/Z lattice R_W (double Delta_fl and double lead delta taken as exact, exact retarded roots).
Does not import wz_rgf, wz_variants, stationary solvers or any mutant."""
import json
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 60
dx = 1 / 8
EP = [4.934802200544679, 19.739208802178716, 44.41321980490211]   # stationary.transverse_basis(1.0, 3).E_perp (doubles)
v0 = 4.0
C2 = dict(mu=(0.5, -0.25, 1, 0), alpha=(1, 2, -1.5, 0.5), beta=(3, -1, 2.5, 1.75), gamma=(1, 0.5, -2, 1.25))
C3 = dict(mu=(0.5, -0.25, 1, 0), a=(1, 2, -1.5, 0.5), a2=(0.25, -0.5, 0.75, 1), b12=(3, -1, 2.5, 1.75), b23=(-2, 1.5, 0.5, -1),
          b13=(1, 0, -1.25, 2), g12=(1, 0.5, -2, 1.25), g23=(-0.5, 1, 0.75, -1.5))
ENER = {"e2a": (2, 12.0), "e2b": (2, 30.0), "e3a": (3, 12.0), "e3b": (3, 30.0), "e3c": (3, 50.0)}


def V(ny, j, cls, s=1.0):
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
        E = lambda a, b: np.eye(3)[:, [a]] @ np.eye(3)[[b], :]
        S = lambda a, b: E(a, b) + E(b, a)
        K = lambda a, b: 1j * (E(b, a) - E(a, b))
        M = c["mu"][k] * np.eye(3) + c["a"][k] * np.diag([1.0, 0, -1]) + c["a2"][k] * np.diag([1.0, -2, 1]) + 0j
        if cls in "RH":
            M = M + c["b12"][k] * S(0, 1) + c["b23"][k] * S(1, 2) + c["b13"][k] * S(0, 2)
        if cls == "H":
            M = M + c["g12"][k] * K(0, 1) + c["g23"][k] * K(1, 2)
    return v0 * s * M          # s and v0 dyadic: exact


def root(d):
    c = 1 - d
    return c + 1j * mp.sqrt(1 - c * c) if d > 0 else c - mp.sqrt(c * c - 1)


def tomp(A):
    return mp.matrix([[mp.mpc(complex(x).real, complex(x).imag) for x in row] for row in A])


def fro(A):
    return mp.sqrt(mp.fsum(abs(A[i, j]) ** 2 for i in range(A.rows) for j in range(A.cols)))


def comm_rel(A, B):
    return fro(A * B - B * A) / (fro(A) * fro(B))


def blk(G, a, b, ny):
    return G[a * ny:(a + 1) * ny, b * ny:(b + 1) * ny]


def cond2(A):
    a = np.array([[complex(A[i, j]) for j in range(A.cols)] for i in range(A.rows)])
    return float(np.linalg.cond(a))


def assemble(Dl, Lam, n, right_lead):
    """M/t of the first n slices: diag 2 Delta_j - 2I, off-diag I, + Lambda on slice 1, + Lambda on slice n if right_lead.
    The two lead terms are added as separate operations (for n = 1 both land on the single slice)."""
    ny = Lam.rows
    M = mp.matrix(n * ny, n * ny)
    for j in range(n):
        Bj = 2 * Dl[j] - 2 * mp.eye(ny)
        for a in range(ny):
            for b in range(ny):
                M[j * ny + a, j * ny + b] = Bj[a, b]
        if j + 1 < n:
            for a in range(ny):
                M[j * ny + a, (j + 1) * ny + a] = 1
                M[(j + 1) * ny + a, j * ny + a] = 1
    for a in range(ny):
        M[a, a] += Lam[a, a]                              # left lead on slice 1
    if right_lead:
        for a in range(ny):
            M[(n - 1) * ny + a, (n - 1) * ny + a] += Lam[a, a]   # right lead on slice n
    return M


def case(eid, cls, N, s=1.0):
    ny, E = ENER[eid]
    Ep = np.array(EP[:ny])
    I = mp.eye(ny)
    Dl_fl = [(E * np.eye(ny) - np.diag(Ep) - V(ny, j, cls, s)) * (dx * dx) for j in range(1, N + 1)]
    Dl = [tomp(D) for D in Dl_fl]
    dlead = [(E - e) * (dx * dx) for e in Ep]
    Lam = mp.diag([root(mp.mpf(float(d))) for d in dlead])
    M = assemble(Dl, Lam, N, True)
    G = M ** -1
    res_full = max(abs(x) for x in (M * G - mp.eye(N * ny)))
    out = dict(case=f"{eid}/{cls}/N{N}" + ("" if s == 1 else f"/s=2^{int(round(np.log2(s)))}"), energy=eid, cls=cls, N=N, s=s,
               n_open=int(sum(d > 0 for d in dlead)), ny=ny, cond_M=cond2(M), norm2_Ghat=float(np.linalg.norm(
                   np.array([[complex(G[i, j]) for j in range(G.cols)] for i in range(G.rows)]), 2)), mp_residual=float(res_full))
    out["kappa_DeltaLambda_max"] = float(max(comm_rel(D, Lam) if fro(D) > 0 else 0 for D in Dl))
    # prefixes j = 1..N-1 (no right lead)
    W, P, Gj1, G1j, G11, condPre, condIW = {}, {}, {}, {}, {}, {}, {}
    for j in range(1, N):
        Mj = assemble(Dl, Lam, j, False)
        Gp = Mj ** -1
        Gjj = blk(Gp, j - 1, j - 1, ny)
        W[j] = I + Gjj ** -1
        P[j] = -Gjj
        Gj1[j], G1j[j], G11[j] = blk(Gp, j - 1, 0, ny), blk(Gp, 0, j - 1, ny), blk(Gp, 0, 0, ny)
        condPre[j] = cond2(Mj)
        condIW[j] = cond2(I - W[j])
    # algebraic closing identity Z = 2 Delta_N + (Lambda - I) + (I - W_{N-1})^{-1} W_{N-1}  vs  Z = Ghat_NN^{-1}
    Zi = blk(G, N - 1, N - 1, ny)
    Z = Zi ** -1
    if N >= 2:
        Zf = 2 * Dl[N - 1] + (Lam - I) + (I - W[N - 1]) ** -1 * W[N - 1]
        out["closing_identity_err"] = float(max(abs(x) for x in (Zf - Z)))
        # prefix chain identities (left/right multiplication, old values) vs direct prefix inverses
        err = 0
        for j in range(2, N):
            err = max(err, max(abs(x) for x in (P[j] * Gj1[j - 1] - Gj1[j])), max(abs(x) for x in (G1j[j - 1] * P[j] - G1j[j])),
                      max(abs(x) for x in (G11[j - 1] - G1j[j - 1] * P[j] * Gj1[j - 1] - G11[j])))
        out["prefix_identity_err"] = float(err)
    else:
        out["N1_identity_err"] = float(max(abs(x) for x in (2 * Dl[0] + 2 * (Lam - I) - Z)))
    out["cond_Z"] = cond2(Z)
    out["cond_prefix_max"] = max(condPre.values()) if condPre else None
    out["cond_I_minus_W_max"] = max(condIW.values()) if condIW else None
    out["normP_max"] = max(cond2(P[j]) * 0 + float(np.linalg.norm(np.array([[complex(P[j][a, b]) for b in range(ny)] for a in range(ny)]), 2)) for j in P) if P else None
    kW = [float(comm_rel(W[j], W[j + 1])) for j in range(1, N - 1)]
    kC = [float(comm_rel(P[j], Gj1[j - 1])) for j in range(2, N)]
    out["kappa_W"] = kW
    out["kappa_chain"] = kC
    out["kappa_Z"] = float(comm_rel(Zi, Gj1[N - 1])) if N >= 2 else None
    # exact detection conditions for M5 (step j >= 3) and M8 (closing), relative to the correct term
    m5 = []
    for j in range(3, N):
        A = G1j[j - 1] * P[j] * Gj1[j - 1]; B = Gj1[j - 1] * P[j] * G1j[j - 1]
        m5.append(float(fro(A - B) / fro(A)))
    out["cond_M5"] = m5
    if N >= 2:
        A = G1j[N - 1] * Zi * Gj1[N - 1]; B = Gj1[N - 1] * Zi * G1j[N - 1]
        out["cond_M8"] = float(fro(A - B) / fro(A))
        out["chain_asym"] = float(fro(Gj1[N - 1] - G1j[N - 1]) / fro(Gj1[N - 1]))
    return out


rows = []
for eid in ENER:
    for cls in "RHD":
        for N in (1, 2, 3, 4):
            rows.append(case(eid, cls, N))
for eid in ("e2b", "e3b"):
    for N in (3, 4):
        for e in (-4, -8, -12):
            rows.append(case(eid, "R", N, 2.0 ** e))
json.dump(rows, open(sys.argv[1], "w"), indent=1)
print(len(rows), "unique configurations evaluated")
