r"""Driver of the W/Z coupling-scaling experiment v2 (PREDICTIONS_coupling_scaling_v2.md, Revision 3.7).

STATUS: implemented, NOT RUN. The main run is gated (`gate`): it refuses to start unless
  * REVIEW_coupling_scaling_v2.md exists, carries an R0 sign-off line, and quotes the manifest.sha256 hashes;
  * results_coupling_scaling_v2/manifest.sha256 matches the committed manifests;
  * previous_files.sha256 matches the current files (§8.3);
  * every open specification gap that the driver cannot resolve without changing the frozen protocol is resolved
    (currently: the §4.7 4(b) "duplicate set" is undefined → `DUPLICATE_SET = None` → the gate refuses).
Tests exercise the pipeline only on fixture configurations (`run_configuration(..., grid=...)`), never on the manifest.

Per (configuration, s) (§4.1, §4.3–§4.7):
  * double inputs exactly as the W/Z code sees them; §4.7 item 1 (Δ_fl off-diagonals bitwise −s·V_off(1)·dx²);
  * WZ solver `wz_rgf.wz_endpoint_blocks`/`wz_fisher_lee`, the instrumented copy (item 6 bit identity), the double
    mutants of `wz_mutants.py`; production `stationary.scattering_matrix` (orig; reported only);
  * validated references (coupling_scaling_v2_refs) of the full and prefix systems on L, R_W and R_orig, with the
    (b) dps-100 and (c) lu_solve checks on their preregistered subsets; INVALID on any failure (item 4);
  * items 2, 3, 5; oracle (§4.2) and transcription (§4.3) on R_W and L; the three §4.3 agreement checks and the
    first-deviation check; η_X, ν_X, τ_abs,m, ν_X,m, SNR*_m, Mag_m;
  * heuristic bounds (τ_end, τ_pre, τ_W, τ_P, τ_S, U/Q/symmetry; inherited from coupled_blocks.bounds/tau_S, §4.6) and
    the validated reference bounds used for "reliably measured" (§6).
Per configuration: §5.2 Taylor coefficients on L (dps 100/200), zero test, consistency checks (i)–(iii) with the
B_tail/B_prec series check and its dps-200 escalation; §5.3 regression; §5.4 slopes, crossover, regimes; §5.5 controls.
Per set: §6 decisions (HO decides; XP is reported only), §4.7 item 7 determinism, §8.1 output files.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import mpmath  # noqa: E402

import coupling_scaling_v2_analysis as AN  # noqa: E402
import coupling_scaling_v2_certify as C  # noqa: E402
import coupling_scaling_v2_refs as R  # noqa: E402
import wz_mutants_mp as TR  # noqa: E402
import wz_mutants_oracle as O  # noqa: E402

OUT_DIR = os.path.join(HERE, "results_coupling_scaling_v2")
REVIEW = os.path.join(HERE, "REVIEW_coupling_scaling_v2.md")
DPS = 60
EPS = 2.0 ** -52
C_SAFE = 10.0
FAIL_THR = 1e-10
JUDGE_THR = 1e-11
DX = 1.0 / 8
T = 32.0
GRID = sorted(C._GRID)                                            # 81 points, increasing
B_SET_BASE = frozenset({Fraction(1, 2**k) for k in range(0, 41, 4)} | {Fraction(3, 2**(k + 2)) for k in range(0, 40, 4)}
                       | {Fraction(1), Fraction(1, 2**40)})       # choice C-5: k ≡ 0 (mod 4) for both grid families
C_SET = frozenset({Fraction(1), Fraction(1, 2**12), Fraction(1, 2**40)})
DUPLICATE_SET = None                                              # §4.7 4(b): undefined in the frozen spec (gap G-1)
BLOCKS = ("G_11", "G_N1", "G_1N", "G_NN")
KEYS = {"G_11": "Ghat11", "G_N1": "GhatN1", "G_1N": "Ghat1N", "G_NN": "GhatNN"}
PURE_ORDER = ("M1", "M2", "M3", "M5", "M6", "M7", "M8", "M15")
STRUCTURAL = ("M4", "M4L", "M4R", "M10", "M13", "M14", "M16")
SERIES_MUTANTS = PURE_ORDER + ("M9", "M12")
EXPECTED_PAIRS = {"HO": {"R": 312, "H": 312}, "XP": {"R": 156, "H": 156}}


class GateRefused(RuntimeError):
    pass


class Invalid(RuntimeError):
    pass


# =========================================================================================================== gate
def _sha256_file(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def gate(out_dir=OUT_DIR, review=REVIEW):
    """Refuse the main run unless R0, the manifest hash, the previous-files hash and the spec gaps are in order."""
    problems = []
    if DUPLICATE_SET is None:
        problems.append("§4.7 4(b) duplicate set is undefined in the frozen spec (gap G-1): needs an R0/spec resolution")
    mh = os.path.join(out_dir, "manifest.sha256")
    hashes = {}
    if not os.path.exists(mh):
        problems.append("manifest.sha256 missing")
    else:
        for line in open(mh):
            h, fn = line.split()
            hashes[fn] = h
            p = os.path.join(out_dir, fn)
            if not os.path.exists(p) or _sha256_file(p) != h:
                problems.append(f"manifest hash mismatch: {fn}")
    pf = os.path.join(out_dir, "previous_files.sha256")
    if not os.path.exists(pf):
        problems.append("previous_files.sha256 missing")
    else:
        for line in open(pf):
            h, rel = line.rstrip("\n").split("  ", 1)
            p = os.path.join(HERE, rel)
            if not os.path.exists(p) or _sha256_file(p) != h:
                problems.append(f"previous file changed or missing: {rel}")
    if not os.path.exists(review):
        problems.append("REVIEW_coupling_scaling_v2.md (R0 sign-off) missing")
    else:
        txt = open(review).read()
        if "R0 SIGN-OFF" not in txt:
            problems.append("REVIEW_coupling_scaling_v2.md has no 'R0 SIGN-OFF' line")
        for fn, h in hashes.items():
            if h not in txt:
                problems.append(f"REVIEW does not quote the hash of {fn}")
    if problems:
        raise GateRefused("; ".join(problems))
    return True


# ===================================================================================================== helpers
def mpf(x):
    x = Fraction(x)
    return mpmath.mpf(x.numerator) / x.denominator


def mpc_pair(z):
    return mpmath.mpc(mpf(z[0]), mpf(z[1]))


def mp_exact(A):
    return mpmath.matrix([[mpc_pair(z) for z in row] for row in A])


def mp_np(A):
    A = np.asarray(A)
    return mpmath.matrix([[mpmath.mpc(complex(x).real, complex(x).imag) for x in row] for row in A])


def np_mp(M):
    return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])


def np_exact(A):
    return np.array([[complex(float(z[0]), float(z[1])) for z in row] for row in A])


def fro(M):
    return mpmath.sqrt(mpmath.fsum(abs(M[i, j]) ** 2 for i in range(M.rows) for j in range(M.cols)))


def maxabs(M):
    return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))


def norm2(M):
    return float(np.linalg.norm(np_mp(M), 2))


def cond2(X):
    return float(np.linalg.cond(X if isinstance(X, np.ndarray) else np_mp(X)))


def blk(G, a, b, ny):
    """1-indexed block."""
    return G[(a - 1) * ny:a * ny, (b - 1) * ny:b * ny]


def endpoint(G, name, N, ny):
    a, b = {"G_11": (1, 1), "G_N1": (N, 1), "G_1N": (1, N), "G_NN": (N, N)}[name]
    return blk(G, a, b, ny)


def cross_fro(A, B):
    """Frobenius distance of two matrices that may live in different mp contexts (entrywise)."""
    return mpmath.sqrt(mpmath.fsum(abs(mpmath.mpc(A[i, j]) - mpmath.mpc(B[i, j])) ** 2
                                   for i in range(B.rows) for j in range(B.cols)))


def _hash_arrays(arrs):
    h = hashlib.sha256()
    for a in arrs:
        h.update(np.ascontiguousarray(np.asarray(a, dtype=complex)).tobytes())
    return h.hexdigest()


# ================================================================================================ configuration
def config_from_manifest(rec):
    V1 = [tuple(tuple((Fraction(z[0]), Fraction(z[1])) for z in row) for row in B) for B in rec["V1"]]
    return dict(set=rec["set"], idx=rec["idx"], ny=rec["N_y"], cls=rec["cls"], N=rec["N"], E=Fraction(rec["E"]),
                n_open=rec["n_open"], V1=V1, E_perp=[float.fromhex(x) for x in rec["E_perp"]],
                cert=rec["certificate"])


# ======================================================================================================= one point
def _solvers():
    import wz_rgf as wz                                            # noqa: E402  (read-only)
    import wz_mutants as wzm                                       # noqa: E402  (read-only)
    from stationary import Device, device_hamiltonian_blocks, scattering_matrix, transverse_basis  # noqa: E402
    import coupled_blocks as CB                                    # noqa: E402  (read-only; inherited τ formulas)
    return wz, wzm, Device, device_hamiltonian_blocks, scattering_matrix, transverse_basis, CB


def double_inputs(cfg, s):
    ny, N, cls = cfg["ny"], cfg["N"], cfg["cls"]
    sd = float(s)
    if Fraction(sd) != s:
        raise Invalid(f"grid point {s} is not exactly representable")
    Epd = np.asarray(cfg["E_perp"], dtype=float)
    V1d = np.array([[[complex(float(z[0]), float(z[1])) for z in r] for r in B] for B in cfg["V1"]])
    if cls in "RD":
        V1d = V1d.real.copy()
    Vs = sd * V1d
    E = float(cfg["E"])
    Delta = (E * np.eye(ny) - np.diag(Epd) - Vs) * (DX * DX)
    dlead = (E - Epd) * (DX * DX)
    return dict(sd=sd, Epd=Epd, V1d=V1d, Vs=Vs, Delta=Delta, dlead=dlead, E=E)


def item1(cfg, s, Delta):
    """§4.7 item 1: off-diagonals bitwise −s·V_off(1)·dx²; returns (ok, max diagonal rounding of the s-part)."""
    ny, N = cfg["ny"], cfg["N"]
    ok, diag = True, Fraction(0)
    Ep = [Fraction(e) for e in cfg["E_perp"]]
    for j in range(N):
        for a in range(ny):
            for b in range(ny):
                z = complex(Delta[j][a, b])
                got = (Fraction(z.real), Fraction(z.imag))
                v = cfg["V1"][j][a][b]
                if a != b:
                    want = (-s * v[0] / 64, -s * v[1] / 64)
                    ok = ok and got == want
                else:
                    exact = (cfg["E"] - Ep[a] - s * v[0]) / 64
                    diag = max(diag, abs(got[0] - exact))
    return ok, diag


def lattices(cfg, s, di, CBmods):
    """Exact matrices and ε_A of the full and prefix systems of L, R_W and R_orig."""
    wz, wzm, Device, dhb, scattering_matrix, tb, CB = CBmods
    ny, N = cfg["ny"], cfg["N"]
    Ep_q = [C.import_binary64(e) for e in cfg["E_perp"]]
    Dq = [C.import_matrix(di["Delta"][j]) for j in range(N)]
    dq = [C.import_binary64(float(x)) for x in di["dlead"]]
    keys = ["full"] + list(range(1, N))
    kind = lambda k: ("full", None) if k == "full" else ("prefix", k)  # noqa: E731
    gamL = [1 - (cfg["E"] - e) * C.DX2 for e in Ep_q]
    gamW = [1 - d for d in dq]
    out = {"L": {}, "R_W": {}, "R_orig": {}}
    for k in keys:
        kd, j = kind(k)
        out["L"][k] = (C.assemble_L(cfg["E"], Ep_q, cfg["V1"], s, kd, j), C.eps_A_L(cfg["E"], Ep_q, N, kd))
        out["R_W"][k] = (C.assemble_R_W(Dq, dq, kd, j), C.eps_A(gamW, N, kd))
    return out, dict(Ep_q=Ep_q, Dq=Dq, dq=dq, gamL=gamL, gamW=gamW)


def orig_lattice(cfg, di, CBmods, out):
    wz, wzm, Device, dhb, scattering_matrix, tb, CB = CBmods
    ny, N = cfg["ny"], cfg["N"]
    dev = Device(xs=DX * np.arange(N), dx=DX, basis=tb(1.0, ny), V=di["Vs"], V_lead_max=0.0)
    rB = scattering_matrix(di["E"], dev, min_closed=0)
    H, tB = dhb(dev)
    if tB != T:
        raise Invalid("production hopping differs from t = 32")
    A_orig = np.array([di["E"] * np.eye(ny, dtype=complex) - H[j] for j in range(N)])
    sig = np.asarray(rB["modes"]["sigma"])
    Aq = [C.import_matrix(A_orig[j]) for j in range(N)]
    sq = [C.import_complex128(complex(x)) for x in sig]
    for k in ["full"] + list(range(1, N)):
        kd, j = ("full", None) if k == "full" else ("prefix", k)
        out["R_orig"][k] = (C.assemble_R_orig(Aq, sq, kd, j), Fraction(0))
    return rB, dict(gamma_code=np.asarray(rB["modes"]["gamma"], dtype=float))


def refs_for(mats, do_b, do_c):
    return {k: R.reference(A, e, do_b, do_c) for k, (A, e) in mats.items()}


def ghat_mp(ref):
    return mp_exact(R.ghat_from_Y(ref["Y"]))


def lam_m1_mp(gammas):
    return [mpc_pair(C.lead_enclosure(g)["lam_mid"]) - 1 for g in gammas]


def lam_m1_exact(gammas):
    return [(C.lead_enclosure(g)["lam_mid"][0] - 1, C.lead_enclosure(g)["lam_mid"][1]) for g in gammas]


def delta_mp_L(cfg, s):
    ny = cfg["ny"]
    Ep = [C.import_binary64(e) for e in cfg["E_perp"]]
    out = []
    for B in cfg["V1"]:
        out.append([[mpmath.mpc(mpf(((cfg["E"] - Ep[a]) if a == b else Fraction(0)) / 64 - s * B[a][b][0] / 64),
                                mpf(-s * B[a][b][1] / 64)) for b in range(ny)] for a in range(ny)])
    return out


def delta_mp_W(Dq):
    return [[[mpc_pair(z) for z in row] for row in B] for B in Dq]


# ---------------------------------------------------------------------------- first-deviation closed forms (§4.3)
FIRST_DEV = {"M1": "G_21^(2)", "M2": "G_12^(2)", "M3": "G_21^(2)", "M4": "G11^(2)", "M4L": "G11^(2)", "M4R": "G11^(2)",
             "M5": "G11^(3)", "M6": "GN1", "M7": "G1N", "M8": "G11", "M9": "P_2", "M10": "Z", "M12": None,
             "M13": "G11^(2)", "M14": "G11^(1)", "M15": "GN1", "M16": "G11"}


def first_dev_closed_form(m, N, ny, G, Gp, lm1):
    """The prefix-level closed form of the first deviating intermediate: the changed line of §4.2 evaluated on the
    prefix references (all its inputs are still exact at that point). Returns (name, value) or (None, None) when the
    change is identical by algebra at this N (M5 at N = 3, M8 at N = 2) or an identity (M12)."""
    P = lambda j: -blk(Gp[j], j, j, ny)                 # noqa: E731
    chi = lambda j: blk(Gp[j], j, 1, ny)                # noqa: E731
    psi = lambda j: blk(Gp[j], 1, j, ny)                # noqa: E731
    G1 = blk(Gp[1], 1, 1, ny) if 1 in Gp else None
    GNN = blk(G, N, N, ny)
    if m in ("M1", "M3"):
        return "G_21^(2)", G1 * P(2)
    if m == "M2":
        return "G_12^(2)", P(2) * G1
    if m == "M4":
        return "G11^(2)", G1 - G1 * P(2) ** 3 * G1
    if m in ("M4L", "M4R"):
        return "G11^(2)", G1 - G1 * P(2) ** 2 * G1
    if m == "M5":
        if N < 4:
            return None, None
        return "G11^(3)", blk(Gp[2], 1, 1, ny) - chi(2) * P(3) * psi(2)
    if m in ("M6", "M15"):
        return "GN1", -(chi(N - 1) * GNN)
    if m == "M7":
        return "G1N", -(GNN * psi(N - 1))
    if m == "M8":
        if N == 2:
            return None, None
        return "G11", blk(Gp[N - 1], 1, 1, ny) + chi(N - 1) * GNN * psi(N - 1)
    if m == "M9":
        return "P_2", -G1
    if m == "M10":
        return "Z", GNN ** -1 - mpmath.diag(lm1)
    if m == "M12":
        return None, None
    if m == "M13":
        return "G11^(2)", G1 + G1 * P(2) * G1
    if m == "M14":
        return "G11^(1)", -G1
    if m == "M16":
        return "G11", blk(Gp[N - 1], 1, 1, ny) + psi(N - 1) * GNN ** 3 * chi(N - 1)
    raise ValueError(m)


def first_deviation(tr_mut, tr_ok, closed, tol):
    """Returns dict(name, ok): the first trace entry deviating by more than tol, which must be the predicted one and
    match the closed form within tol; with no closed form (identity / identical by algebra) no entry may deviate."""
    ok_map = dict(tr_ok)
    first = None
    for nm, X in tr_mut:
        Y = ok_map.get(nm)
        if Y is None or (X.rows, X.cols) != (Y.rows, Y.cols) or cross_fro(X, Y) > tol:
            first = (nm, X)
            break
    name, val = closed
    if name is None:
        return dict(expected=None, first=None if first is None else first[0], ok=first is None, value=None)
    if first is None or first[0] != name:
        return dict(expected=name, first=None if first is None else first[0], ok=False, value=None)
    v = cross_fro(first[1], val)
    return dict(expected=name, first=name, ok=bool(v <= tol), value=float(v))


# ------------------------------------------------------------------------------------------- the per-point pipeline
def point(cfg, s, mods):
    wz, wzm, Device, dhb, scattering_matrix, tb, CB = mods
    ny, N, cls = cfg["ny"], cfg["N"], cfg["cls"]
    n_tot = N * ny
    inval, anomalies = [], []
    with mpmath.workdps(DPS):
        di = double_inputs(cfg, s)
        ok1, diag_round = item1(cfg, s, di["Delta"])
        if not ok1:
            inval.append("§4.7 item 1: off-diagonal Delta_fl not bitwise -s V_off dx^2")
        lead = wz.lead_from_delta(di["dlead"])
        GA = wz.wz_endpoint_blocks(di["Delta"], lead["lam_minus_1"])
        SA = wz.wz_fisher_lee(GA, lead["sin_theta"], lead["open"])
        inst = wzm.wz_instrumented(di["Delta"], lead["lam_minus_1"])
        faithful = all(np.array_equal(GA[k], inst[k]) for k in KEYS.values())
        if not faithful:
            inval.append("§4.7 item 6: instrumented copy not bit-identical")
        mats, aux = lattices(cfg, s, di, mods)
        rB, odat = orig_lattice(cfg, di, mods, mats)
        do_b = s in B_SET_BASE or (DUPLICATE_SET is not None and s in DUPLICATE_SET)
        do_c = s in C_SET
        refs = {lat: refs_for(mats[lat], do_b, do_c) for lat in ("L", "R_W", "R_orig")}
        for lat, rr in refs.items():
            for k, r in rr.items():
                if not r["valid"]:
                    inval.append(f"§4.7 item 4: reference {lat}/{k} not validated ({r['status']}, checks {r.get('checks')})")
        if inval and any("item 4" in x for x in inval):
            return dict(s=str(s), invalid=inval, aborted=True)
        G = {lat: ghat_mp(refs[lat]["full"]) for lat in refs}
        Gp = {lat: {j: ghat_mp(refs[lat][j]) for j in range(1, N)} for lat in refs}
        normG = {lat: norm2(G[lat]) for lat in G}
        errG = {lat: R.err_ghat(refs[lat]["full"]) for lat in refs}
        condM = {lat: cond2(np_exact(mats[lat]["full"][0])) for lat in mats}
        # item 2
        rep = {}
        for lat in ("R_W", "R_orig"):
            rep[lat] = float(fro(G[lat] - G["L"]) / normG["L"])
            if not rep[lat] <= 10 * condM["L"] * EPS:
                inval.append(f"§4.7 item 2: e_glob(Ref({lat}) - Ref(L)) = {rep[lat]:.3e} > 10 cond eps")
        # item 5 (R_W and L)
        idv = {}
        for lat in ("R_W", "L"):
            idv[lat] = float(identities(G[lat], Gp[lat], mats[lat], lat, aux, N, ny) / normG[lat])
            if not idv[lat] <= 1e-50:
                inval.append(f"§4.7 item 5: identities on {lat} violated ({idv[lat]:.2e} ||G||)")
        # W_j / P_j validated references and bounds (R_W)
        W_ref = {j: mp_exact(w) for j, w in R.w_reference(mats["R_W"]["full"][0], refs["R_W"], N, ny).items()}
        wp_bnd = R.w_p_bounds({j: refs["R_W"][j] for j in range(1, N)}, mats["R_W"][1][1] if N > 1 else 0, N)
        bstruct = dict(cond_M=condM["R_W"], Z=blk(G["R_W"], N, N, ny) ** -1, norm2_G=normG["R_W"],
                       prefix={j: dict(W=W_ref[j], P=-blk(Gp["R_W"][j], j, j, ny),
                                       cond=cond2(np_exact(mats["R_W"][j][0]))) for j in range(1, N)})
        n_open = int(np.sum(lead["open"]))
        bW = CB.bounds(bstruct, N, ny, n_open)
        A = solver_A(cfg, GA, SA, inst, G["R_W"], Gp["R_W"], W_ref, wp_bnd, refs["R_W"], aux["gamW"], bW, normG["R_W"],
                     errG["R_W"], mods, N, ny, lead)
        B = solver_B(cfg, rB, G["R_orig"], normG["R_orig"], odat, mods, N, ny)
        # mutants
        muts = mutants(cfg, s, di, lead, inst, G, Gp, aux, bW, normG, condM, mods, N, ny, inval)
    rec = dict(s=str(s), invalid=inval, aborted=False, faithful=faithful, item1_ok=ok1, diag_rounding=float(diag_round),
               rep_diff=rep, rep_bound=10 * condM["L"] * EPS, identities=idv, cond=condM, normG=normG,
               ref_meta={lat: {str(k): dict(dps=r["dps"], rel=float(r["rel"]), escalated=r["escalated"],
                                            checks={c: v for c, v in r["checks"].items() if c in ("b", "c")})
                               for k, r in rr.items()} for lat, rr in refs.items()},
               tau_end=bW["tau_end"], K_Schur=bW["K_Schur"], g=bW["g"], judgeable=bW["tau_end"] <= JUDGE_THR,
               A=A, B=B, mutants=muts,
               hash=_hash_arrays([GA[k] for k in KEYS.values()] + [SA["S"]] + [rB["S"]]
                                 + [m["double_hash_arrays"] for m in muts.values() if "double_hash_arrays" in m]))
    for m in muts.values():
        m.pop("double_hash_arrays", None)
    return rec


def identities(G, Gp, mats, lat, aux, N, ny):
    """§4.7 item 5: prefix chain identities and the closing identities (max entry, Ĝ units)."""
    I = mpmath.eye(ny)
    P = lambda j: -blk(Gp[j], j, j, ny)                 # noqa: E731
    chi = lambda j: blk(Gp[j], j, 1, ny)                # noqa: E731
    psi = lambda j: blk(Gp[j], 1, j, ny)                # noqa: E731
    G11p = lambda j: blk(Gp[j], 1, 1, ny)               # noqa: E731
    e = mpmath.mpf(0)
    for j in range(2, N):
        e = max(e, maxabs(P(j) * chi(j - 1) - chi(j)), maxabs(psi(j - 1) * P(j) - psi(j)),
                maxabs(G11p(j - 1) - psi(j - 1) * P(j) * chi(j - 1) - G11p(j)))
    GNN = blk(G, N, N, ny)
    e = max(e, maxabs(G11p(N - 1) + psi(N - 1) * GNN * chi(N - 1) - blk(G, 1, 1, ny)),
            maxabs(-(GNN * chi(N - 1)) - blk(G, N, 1, ny)), maxabs(-(psi(N - 1) * GNN) - blk(G, 1, N, ny)))
    return e


def solver_A(cfg, GA, SA, inst, G, Gp, W_ref, wp_bnd, refs, gam, bW, normG, errG, mods, N, ny, lead):
    CB = mods[-1]
    cls = cfg["cls"]
    sq = math.sqrt(ny)
    rb_end = float(sq * errG) / (normG - float(errG)) if errG is not None else None
    blocks = {}
    for name in BLOCKS:
        v = CB.block_err(GA[KEYS[name]], endpoint(G, name, N, ny), normG)
        v["ref_bound"] = rb_end
        v["e_rel_bound"] = bW["tau_end"] * normG / v["ref_fro"] if v["ref_fro"] > 0 else None
        blocks[name] = v
    srec = R.s_reference(refs["full"], gam, N, ny)
    no = len(srec["open"])
    Sd = SA["S"]
    dS, dS_bound = 0.0, 0.0
    Sref = mpmath.matrix(2 * no, 2 * no)
    for (a, b), z in srec["S"].items():
        Sref[a, b] = mpc_pair(z)
        dS = max(dS, abs(complex(Sd[a, b]) - complex(Sref[a, b])))
        bd = srec["bound"][(a, b)]
        dS_bound = None if (bd is None or dS_bound is None) else max(dS_bound, float(bd))
    Ur, Qr = CB.UQ_mp(Sref)
    U, Q = CB.UQ_np(Sd)
    smax = float(max(lead["sin_theta"][n] for n in srec["open"])) if no else 0.0
    tS = CB.tau_S(dict(norm2_G=normG), bW, smax)
    res = dict(blocks=blocks, dS_max=dS, dS_ref_bound=dS_bound, U=U, Q=Q, U_ref=float(Ur), Q_ref=float(Qr),
               tau_end=bW["tau_end"], tau_S=tS, U_bound=float(Ur) + 2 * (2 * no) * tS, n_open=no,
               sym_residual=float(np.max(np.abs(GA["Ghat1N"] - GA["GhatN1"].T))), sym_bound=2 * bW["tau_end"] * normG,
               sym_applicable=cls in "RD")
    res["Q_metric"] = abs(Q - float(Qr))
    res["Q_bound"] = (2 * tS) if cls == "H" else float(Qr) + 2 * tS
    res["Q_value_judged"] = res["Q_metric"] if cls == "H" else Q
    inter, tr = [], dict(inst["trace"])
    for nm, X in inst["trace"]:
        if nm.startswith("W_") or nm.startswith("P_"):
            j = int(nm[2:])
            ref = W_ref[j] if nm[0] == "W" else -blk(Gp[j], j, j, ny)
            e = wp_bnd[j]["err_W" if nm[0] == "W" else "err_P"]
            inter.append(dict(name=nm, kind=nm[0], value=float(fro(mp_np(X) - ref)),
                              bound=bW["tau_W"][j] if nm[0] == "W" else bW["tau_P"][j],
                              ref_bound=None if e is None else float(sq * e)))
        elif nm.startswith("G11^(") or nm.startswith("G_"):
            j = int(nm.split("^(")[1][:-1])
            if nm.startswith("G11^"):
                ref = blk(Gp[j], 1, 1, ny)
            elif nm.startswith(f"G_{j}1"):
                ref = blk(Gp[j], j, 1, ny)
            else:
                ref = blk(Gp[j], 1, j, ny)
            nj = norm2(Gp[j])
            e = R.err_ghat(refs[j])
            inter.append(dict(name=nm, kind="prefix", value=float(fro(mp_np(X) - ref) / nj), bound=bW["tau_pre"][j],
                              ref_bound=None if e is None else float(sq * e) / (nj - float(e))))
    resid = [dict(name="W_1 - I", rho=CB.rho(tr["W_1"] - np.eye(ny), tr["G11^(1)"]))]
    for j in range(2, N):
        resid.append(dict(name=f"I - W_{j}", rho=CB.rho(np.eye(ny) - tr[f"W_{j}"], tr[f"P_{j}"])))
    resid.append(dict(name="Z", rho=CB.rho(tr["Z"], tr["Zi"])))
    res.update(intermediates=inter, residuals=resid)
    return res


def solver_B(cfg, rB, G, normG, odat, mods, N, ny):
    """Production orig solver vs Ref(R_orig): accuracy reported only (no decision weight)."""
    CB = mods[-1]
    GB = {"G_11": T * rB["G"]["G11"], "G_N1": T * rB["G"]["GN1"], "G_1N": T * rB["G"]["G1N"], "G_NN": T * rB["G"]["GNN"]}
    blocks = {name: CB.block_err(GB[name], endpoint(G, name, N, ny), normG) for name in BLOCKS}
    g = odat["gamma_code"]
    o = [n for n in range(ny) if g[n] > 0]
    f = [mpmath.sqrt(mpmath.mpf(float(g[n]))) for n in o]
    no = len(o)
    S = mpmath.matrix(2 * no, 2 * no)
    pos = {(0, 0): (1, 1), (0, 1): (1, N), (1, 0): (N, 1), (1, 1): (N, N)}
    for (sa, sb), (ba, bb) in pos.items():
        X = blk(G, ba, bb, ny)
        for i in range(no):
            for k in range(no):
                S[sa * no + i, sb * no + k] = 1j * f[i] * f[k] / T * X[o[i], o[k]] - (1 if (sa == sb and i == k) else 0)
    dS = float(maxabs(mp_np(rB["S"]) - S))
    return dict(blocks=blocks, max_e_glob=max(v["e_glob"] for v in blocks.values()), dS_max=dS)


def mutants(cfg, s, di, lead, inst, G, Gp, aux, bW, normG, condM, mods, N, ny, inval):
    wz, wzm = mods[0], mods[1]
    cls = cfg["cls"]
    n_tot = N * ny
    be = O.MpBackend(mpmath.mp)
    lm1 = {"R_W": lam_m1_mp(aux["gamW"]), "L": lam_m1_mp(aux["gamL"])}
    dl = {"R_W": delta_mp_W(aux["Dq"]), "L": delta_mp_L(cfg, s)}
    tr_ok = {lat: TR.mp_endpoint_blocks(dl[lat], lm1[lat], None, dps=DPS) for lat in ("R_W", "L")}
    tol = {lat: mpmath.mpf(10) ** -50 * normG[lat] for lat in ("R_W", "L")}
    for lat in ("R_W", "L"):
        for name in BLOCKS:
            d = cross_fro(tr_ok[lat][KEYS[name]], endpoint(G[lat], name, N, ny))
            if d > tol[lat]:
                inval.append(f"§4.3 transcription correct code != Ref({lat}) at {name} ({float(d):.2e})")
    tau_abs = bW["tau_end"] * normG["R_W"]
    Kmax_c = max(condM["R_W"], bW["K_Schur"])
    out = {}
    for m in wzm.MUTANTS:
        if not O.on_path(m, N) or m == "M11":
            out[m] = dict(mutant=m, status="off_path")
            continue
        r = wzm.wz_instrumented(di["Delta"], lead["lam_minus_1"], mutant=m)
        des = O.DESIGNATED[m]
        rec = dict(mutant=m, status="computed", designated=des)
        orc = {lat: O.oracle(m, N, ny, G[lat], Gp[lat], lm1[lat], be) for lat in ("R_W", "L")}
        trm = {lat: TR.mp_endpoint_blocks(dl[lat], lm1[lat], m, dps=DPS) for lat in ("R_W", "L")}
        Xref = {lat: endpoint(G[lat], des, N, ny) for lat in ("R_W", "L")}
        Xd = mp_np(r[KEYS[des]])
        D = float(fro(Xd - Xref["R_W"]))
        Dst = {lat: float(fro(orc[lat]["X"] - Xref[lat])) for lat in ("R_W", "L")}
        eta = float(fro(mp_np(inst[KEYS[des]]) - Xref["R_W"]))
        fXref = float(fro(Xref["R_W"]))
        nu = max(eta, tau_abs, EPS * fXref)
        Xor = orc["R_W"]["X"]
        mag = float(O.magnitude(m, N, ny, G["R_W"], Gp["R_W"], be, norm2, X_oracle=Xor if m == "M10" else None))
        tau_abs_m = max(tau_abs, C_SAFE * n_tot * EPS * Kmax_c * mag)
        nu_m = max(eta, tau_abs_m, EPS * fXref)
        snr = Dst["R_W"] / nu if nu > 0 else None
        # agreement checks (§4.3)
        a_ot = {lat: float(cross_fro(trm[lat][KEYS[des]], orc[lat]["X"])) for lat in ("R_W", "L")}
        for lat in ("R_W", "L"):
            for k, v in orc[lat]["secondary"].items():
                a_ot[lat] = max(a_ot[lat], float(cross_fro(trm[lat][KEYS[k]], v)))
            if not a_ot[lat] <= 1e-50 * normG[lat]:
                inval.append(f"§4.3 oracle vs transcription ({lat}) {m}: {a_ot[lat]:.2e} > 1e-50 ||G||")
        a_do = float(fro(Xd - Xor))
        do_applicable = snr is not None and snr >= 10
        a_dt = max(float(cross_fro(trm["R_W"][KEYS[b]], mp_np(r[KEYS[b]]))) for b in BLOCKS)
        fd = first_deviation(trm["R_W"]["trace"], tr_ok["R_W"]["trace"],
                             first_dev_closed_form(m, N, ny, G["R_W"], Gp["R_W"], lm1["R_W"]), tol["R_W"])
        if not fd["ok"]:
            inval.append(f"§4.3 first-deviation check failed for {m}: {fd}")
        # item 3
        i3 = None
        if Dst["L"] >= 1e-30:
            i3 = abs(Dst["R_W"] - Dst["L"]) / Dst["L"]
            if not i3 <= 1e-6:
                inval.append(f"§4.7 item 3: {m} |D*_RW - D*_L|/D*_L = {i3:.2e} > 1e-6")
        four = {}
        for b in BLOCKS:
            Xr = endpoint(G["R_W"], b, N, ny)
            four[b] = dict(D=float(fro(mp_np(r[KEYS[b]]) - Xr)),
                           Dtr_RW=float(cross_fro(trm["R_W"][KEYS[b]], Xr)),
                           Dtr_L=float(cross_fro(trm["L"][KEYS[b]], endpoint(G["L"], b, N, ny))),
                           Dor_RW=(Dst["R_W"] if b == des else (float(fro(orc["R_W"]["secondary"][b] - Xr))
                                                                 if b in orc["R_W"]["secondary"] else None)))
        rec.update(D=D, Dstar_RW=Dst["R_W"], Dstar_L=Dst["L"], Dstar_tr_RW=float(cross_fro(trm["R_W"][KEYS[des]], Xref["R_W"])),
                   eta=eta, nu=nu, tau_abs=tau_abs, tau_abs_m=tau_abs_m, nu_m=nu_m, Mag=mag, snr=snr, Xref_fro=fXref,
                   agree_oracle_transcription=a_ot, agree_double_oracle=dict(value=a_do, applicable=do_applicable,
                                                                             ok=(not do_applicable) or a_do <= 10 * nu_m),
                   agree_double_transcription=dict(value=a_dt, ok=a_dt <= 10 * nu_m), first_deviation=fd,
                   item3=i3, blocks=four, detected=(D >= 10 * nu) and ((not do_applicable) or a_do <= 10 * nu_m),
                   double_hash_arrays=np.concatenate([r[k].ravel() for k in KEYS.values()]))
        out[m] = rec
    # correct-code noise per block (§5.5)
    out["_noise"] = {b: dict(eta=float(fro(mp_np(inst[KEYS[b]]) - endpoint(G["R_W"], b, N, ny))),
                             eps_ref=EPS * float(fro(endpoint(G["R_W"], b, N, ny)))) for b in BLOCKS}
    out["_noise"]["double_hash_arrays"] = np.zeros(1)
    return out


# ================================================================================================= §5.2 per pair
def L_ctx(cfg, s, kind, j, ctx):
    """L(s) as a ctx matrix with the EXACT lead roots evaluated at ctx precision (§5.2 works on the exact lattice L; the
    rational enclosure midpoint of the certificate is accurate only to 2^-256, which would leave e.g. M9's C_0 at
    ~1e-77 instead of 0, Fact 6)."""
    Ep_q = [C.import_binary64(e) for e in cfg["E_perp"]]
    A = C.assemble_L(cfg["E"], Ep_q, cfg["V1"], s, kind, j)
    q = lambda x: ctx.mpf(x.numerator) / x.denominator                                   # noqa: E731
    M = ctx.matrix([[ctx.mpc(q(z[0]), q(z[1])) for z in r] for r in A])
    ny, N = len(Ep_q), cfg["N"]
    n = N if kind == "full" else j
    for p_, e in enumerate(Ep_q):
        g = 1 - (cfg["E"] - e) * C.DX2
        enc = C.lead_enclosure(g)
        gm = q(g)
        lam = ctx.mpc(gm, ctx.sqrt(1 - gm * gm)) if enc["branch"] == "open" else ctx.mpc(gm - ctx.sqrt(gm * gm - 1), 0)
        corr = 32 * (lam - ctx.mpc(q(enc["lam_mid"][0]), q(enc["lam_mid"][1])))
        M[p_, p_] += corr
        if kind == "full":
            k = (n - 1) * ny + p_
            M[k, k] += corr
    return M


def series_for(cfg, dps):
    """Taylor coefficients (dps) of Ĝ_full and every prefix on L: {'full': [C_k], j: [C_k]}."""
    ctx = mpmath.MPContext()
    ctx.dps = dps
    N = cfg["N"]
    out = {}
    for k in ["full"] + list(range(1, N)):
        kd, j = ("full", None) if k == "full" else ("prefix", k)
        A0 = L_ctx(cfg, Fraction(0), kd, j, ctx)
        A1 = C.assemble_A1(cfg["V1"], kd, j)
        A1m = ctx.matrix([[ctx.mpc(ctx.mpf(z[0].numerator) / z[0].denominator, ctx.mpf(z[1].numerator) / z[1].denominator)
                           for z in r] for r in A1])
        out[k] = AN.neumann_series(ctx, A0, A1m, 32)
    return ctx, out


def direct_at(cfg, s, dps):
    ctx = mpmath.MPContext()
    ctx.dps = dps
    N = cfg["N"]
    out = {}
    for k in ["full"] + list(range(1, N)):
        kd, j = ("full", None) if k == "full" else ("prefix", k)
        out[k] = 32 * ctx.inverse(L_ctx(cfg, s, kd, j, ctx))
    return ctx, out


def _fro_ctx(ctx, X):
    return ctx.sqrt(ctx.fsum(abs(X[i, j]) ** 2 for i in range(X.rows) for j in range(X.cols)))


def pair_series(cfg, m, ser, normG0):
    """§5.2 zero test and order for pair (cfg, m) from precomputed series at dps 100 and 200."""
    N, ny = cfg["N"], cfg["ny"]
    zs, Cn = {}, {}
    for dps, (ctx, S) in ser.items():
        be = AN.SeriesBackend(ctx)
        X = O.oracle(m, N, ny, S["full"], {j: S[j] for j in range(1, N)}, None, be)["X"]
        Xr = O.reference_block(O.DESIGNATED[m], N, ny, S["full"], be)
        Ck = [x - y for x, y in zip(X, Xr)]
        Cn[dps] = [_fro_ctx(ctx, c) for c in Ck]
        zs[dps] = [c / normG0 for c in Cn[dps]]
    with mpmath.workdps(220):                     # the zero test compares at 1e-20 relative and 1e-170 absolute
        tests = [AN.zero_test(mpmath.mpf(z1), mpmath.mpf(z2)) for z1, z2 in zip(zs[100], zs[200])]
    return dict(tests=tests, order=AN.identify_order(tests), Cnorm100=[float(x) for x in Cn[100]],
                Cnorm200=[float(x) for x in Cn[200]], Cnorm200_mp=Cn[200])


def series_check(cfg, m, ser, r0, a0s, beta):
    """§5.2 (ii) series vs direct at s*, with B_tail (validated) + B_prec; dps-200 escalation."""
    N, ny = cfg["N"], cfg["ny"]
    mono = O.monomials(m, N)
    Tm = sum(c for c, _ in mono)
    deg = max(d for _, d in mono)
    alpha = AN.alpha_from(a0s)
    B, s_used, x, zeta = AN.b_tail(Tm, deg, alpha, beta, AN.s_star(r0))
    out = dict(T_m=Tm, deg_m=deg, alpha=str(alpha), s_star=str(s_used), x=str(x), B_tail=str(B), levels=[])
    for dps in (100, 200):
        ctx, S = ser[dps]
        dctx, Dd = direct_at(cfg, s_used, dps)
        beS, beD = AN.SeriesBackend(ctx), O.MpBackend(dctx)
        X = O.oracle(m, N, ny, S["full"], {j: S[j] for j in range(1, N)}, None, beS)["X"]
        Xr = O.reference_block(O.DESIGNATED[m], N, ny, S["full"], beS)
        sv = ctx.mpf(s_used.numerator) / s_used.denominator
        SK = None
        for k, (a, b) in enumerate(zip(X, Xr)):
            term = (a - b) * sv ** k
            SK = term if SK is None else SK + term
        Fd = O.oracle(m, N, ny, Dd["full"], {j: Dd[j] for j in range(1, N)}, None, beD)["X"] - \
            O.reference_block(O.DESIGNATED[m], N, ny, Dd["full"], beD)
        diff = ctx.sqrt(ctx.fsum(abs(ctx.mpc(SK[i, j]) - ctx.mpc(Fd[i, j])) ** 2 for i in range(SK.rows) for j in range(SK.cols)))
        mag = O.magnitude(m, N, ny, Dd["full"], {j: Dd[j] for j in range(1, N)}, beD,
                          lambda M_: float(np.linalg.norm(np_mp(M_), 2)))
        with mpmath.workdps(dps):
            Bp = AN.b_prec(dps, mpmath.mpf(mag))
            lhs = mpmath.mpf(diff)
            ok = lhs <= mpf(B) + Bp
        out["levels"].append(dict(dps=dps, diff=float(lhs), B_prec=float(Bp), ok=bool(ok)))
        if ok:
            out["status"] = "pass" if dps == 100 else "pass_after_escalation"
            return out
    out["status"] = "unresolved series inconsistency"
    return out


def a0_and_beta(cert):
    """a_A,0 = 1/(2 r_0,A b_A) from the manifest (b > 0); β from the manifest."""
    a0s = []
    for name, r in cert["matrices"].items():
        b = Fraction(r["b"])
        z = next(p for p in r["points"] if p["s"] == "0/1")
        if b > 0:
            a0s.append(1 / (2 * Fraction(z["r"]) * b))
    return a0s, (Fraction(cert["beta"]) if cert["beta"] is not None else None)


def category(m, cfg):
    N, cls = cfg["N"], cfg["cls"]
    if not O.on_path(m, N) or m == "M11":
        return "off_path"
    if (m == "M5" and N == 3) or (m == "M8" and N == 2):
        return "identical_by_algebra"
    if m == "M12":
        return "equivalent"
    if m in PURE_ORDER and cls == "D":
        return "class_D_order"
    if m == "M9":
        return "perturbation-control"
    if m in STRUCTURAL:
        return "structural"
    return "candidate"


# ============================================================================================= one configuration
def run_configuration(cfg, grid=GRID, mods=None):
    mods = mods or _solvers()
    pts = []
    for s in grid:
        pts.append(point(cfg, s, mods))
        if pts[-1].get("aborted"):
            break
    invalid = [x for p in pts for x in p["invalid"]]
    res = dict(config={k: cfg[k] for k in ("set", "idx", "ny", "cls", "N", "n_open")} | dict(E=str(cfg["E"])),
               points=pts, invalid=invalid, pairs={}, noise_trend={})
    if any(p.get("aborted") for p in pts):
        return res
    svals = [Fraction(p["s"]) for p in pts]
    for b in BLOCKS:
        res["noise_trend"][b] = AN.noise_trend(svals, [p["mutants"]["_noise"][b]["eta"] for p in pts],
                                               [p["mutants"]["_noise"][b]["eps_ref"] for p in pts])
    r0 = cfg["cert"]["r0"]
    a0s, beta = a0_and_beta(cfg["cert"])
    ser = {dps: series_for(cfg, dps) for dps in (100, 200)}
    normG0 = mpmath.mpf(float(np.linalg.norm(np_mp(ser[100][1]["full"][0]), 2)))
    for m in O.MUTANTS:
        cat = category(m, cfg)
        pr = dict(mutant=m, category=cat)
        if cat == "off_path":
            res["pairs"][m] = pr
            continue
        D = [p["mutants"][m]["D"] for p in pts]
        Ds = [p["mutants"][m]["Dstar_RW"] for p in pts]
        snr = [p["mutants"][m]["snr"] for p in pts]
        tid = [10.0 ** -(DPS - 10) * p["normG"]["R_W"] for p in pts]
        pr["all_Dstar_le_tau_id"] = all(d <= t for d, t in zip(Ds, tid))
        pr["apparatus_anomaly"] = any(not (p["mutants"][m]["agree_double_oracle"]["ok"]
                                           and p["mutants"][m]["agree_double_transcription"]["ok"]) for p in pts)
        if m in SERIES_MUTANTS:
            ps = pair_series(cfg, m, ser, normG0)
            pr.update(tests=ps["tests"], Cnorm100=ps["Cnorm100"], Cnorm200=ps["Cnorm200"], order=list(ps["order"]))
            kind, q = ps["order"]
            if m in PURE_ORDER + ("M9",) and ps["tests"][0] == "nonzero":
                res["invalid"].append(f"§5.2 (i): {m} C_0 nonzero (theorem requires zero)")
            if (cat in ("class_D_order", "equivalent")) and not (kind == "all_zero" and pr["all_Dstar_le_tau_id"]):
                res["invalid"].append(f"§5.2 (iii): {m} ({cat}) coefficients/grid discrepancies not zero within tau_id")
            pr["c0_unresolved"] = ps["tests"][0] == "unresolved"
            if cat == "candidate":
                if kind == "q" and q >= 1:
                    pr["category"] = "ordering-meaningful"
                elif kind == "all_zero" and pr["all_Dstar_le_tau_id"]:
                    pr["category"] = "ordering-zero"
                else:
                    pr["category"] = "order-unresolved"
            if kind == "q" and q >= 1:
                pr["q"] = q
                if q + 1 <= AN.K_MAX and ps["Cnorm200"][q] > 0:
                    pr["mu_est"] = ps["Cnorm200"][q + 1] / ps["Cnorm200"][q]
                if pr["category"] in ("ordering-meaningful", "perturbation-control"):
                    pr["series_check"] = series_check(cfg, m, ser, r0, a0s, beta)
                    pr.update(analysis(svals, D, Ds, snr, q, ps["Cnorm200"][q]))
        if cat == "structural":
            pr["hyp4_structural"] = all(p["mutants"][m]["detected"] for p in pts
                                        if p["mutants"][m]["snr"] is not None and p["mutants"][m]["snr"] >= 10)
        if cat in ("class_D_order", "equivalent"):
            pr["hyp4_control"] = all(p["mutants"][m]["D"] <= 10 * p["mutants"][m]["nu"] for p in pts)
        res["pairs"][m] = pr
    return res


def analysis(svals, D, Ds, snr, q, Cq):
    pts = list(zip(svals, D, Ds, snr))
    reg = AN.regression(pts, q)
    phi_e = AN.local_slopes(svals, Ds)
    phi_d = AN.local_slopes(svals, D)
    win = {p[0] for p in AN.window(pts)}
    adequate = AN.window_adequacy(reg, [v for s, v in phi_e.items() if s in win], q)
    # Hyp1 wherever both neighbours have SNR* >= 1e3
    hyp1 = True
    for k in range(1, len(svals) - 1):
        s = svals[k]
        if s in phi_e and s in phi_d and snr[k - 1] is not None and snr[k + 1] is not None \
                and snr[k - 1] >= 1e3 and snr[k + 1] >= 1e3:
            hyp1 = hyp1 and abs(phi_d[s] - phi_e[s]) <= 0.01
    if Fraction(1) in phi_e and Fraction(1) in phi_d and snr[-1] is not None and snr[-1] >= 1e3 \
            and snr[svals.index(Fraction(1, 2))] is not None and snr[svals.index(Fraction(1, 2))] >= 1e3:
        hyp1 = hyp1 and abs(phi_d[Fraction(1)] - phi_e[Fraction(1)]) <= 0.01
    sx, sx2 = AN.crossover(svals, Ds, Cq, q, phi_e)
    return dict(regression=reg, window_adequate=adequate, hyp1=hyp1, hyp2=reg["hyp2"], s_cross=str(sx) if sx else None,
                s_cross_secondary=str(sx2) if sx2 else None, regime=AN.regime_models(svals, Ds, q),
                phi1_exact=phi_e.get(Fraction(1)), phi1_double=phi_d.get(Fraction(1)))


# ============================================================================================== §6 decisions
def decide(set_name, results, determinism_ok=True):
    """§6 decisions per scope (physical: class R + class D controls; algebraic: class H). Returns dict."""
    inv = [x for r in results for x in r["invalid"]]
    if not determinism_ok:
        inv.append("§4.7 item 7: determinism mismatch")
    scopes = {"physical": ("R", "D"), "algebraic": ("H",)}
    out = dict(set=set_name, decides=set_name == "HO", invalid=inv, scopes={})
    for sc, classes in scopes.items():
        fw, fd, di = [], [], []
        rs = [r for r in results if r["config"]["cls"] in classes]
        for r in rs:
            key = f"{set_name}#{r['config']['idx']}"
            for p in r["points"]:
                if p.get("aborted"):
                    continue
                n_fw0 = len(fw)
                A = p["A"]
                judge = p["judgeable"]
                quantities = [(f"endpoint {b}", v["e_glob"], v["ref_bound"], A["tau_end"]) for b, v in A["blocks"].items()]
                quantities += [("dS", A["dS_max"], A["dS_ref_bound"], A["tau_S"])]
                quantities += [(it["name"], it["value"], it["ref_bound"], it["bound"]) for it in A["intermediates"]]
                for nm, val, rb, hb in quantities:
                    if rb is not None and val - rb > FAIL_THR:
                        fw.append(f"{key} s={p['s']}: {nm} {val:.2e} reliably > 1e-10")
                    elif val > FAIL_THR and rb is None:
                        di.append(f"{key} s={p['s']}: {nm} {val:.2e} > 1e-10 not reliably measured")
                    elif judge and hb is not None and val > hb and val <= FAIL_THR:
                        di.append(f"{key} s={p['s']}: anomaly {nm} {val:.2e} > heuristic {hb:.2e}")
                if A["U"] > A["U_bound"] and judge:
                    di.append(f"{key} s={p['s']}: anomaly U")
                if A["Q_value_judged"] > A["Q_bound"] and judge:
                    di.append(f"{key} s={p['s']}: anomaly Q")
                if A["sym_applicable"] and A["sym_residual"] > A["sym_bound"] and judge:
                    di.append(f"{key} s={p['s']}: anomaly symmetry")
                if not judge and len(fw) == n_fw0:
                    di.append(f"{key} s={p['s']}: non-judgeable (tau_end {p['tau_end']:.2e})")
            for b, nt in r["noise_trend"].items():
                if not nt["ok"]:
                    di.append(f"{key}: noise-trend control failed for {b} (slope {nt['slope']})")
            for m, pr in r["pairs"].items():
                c = pr["category"]
                if pr.get("apparatus_anomaly"):
                    di.append(f"{key} {m}: mutant-apparatus anomaly")
                if pr.get("c0_unresolved") and c in ("ordering-meaningful", "order-unresolved", "perturbation-control"):
                    di.append(f"{key} {m}: unresolved C_0")
                if c == "order-unresolved":
                    di.append(f"{key} {m}: order-unresolved")
                if c == "ordering-meaningful":
                    reg = pr.get("regression", {})
                    if not (reg.get("valid") and pr.get("window_adequate")):
                        di.append(f"{key} {m}: not regression-valid or not window-adequate")
                    elif pr.get("hyp2") is False:
                        fd.append(f"{key} {m}: Hyp2 violated")
                    if pr.get("hyp1") is False:
                        fd.append(f"{key} {m}: Hyp1 violated")
                    if pr.get("series_check", {}).get("status") == "unresolved series inconsistency":
                        di.append(f"{key} {m}: unresolved series inconsistency")
                if pr.get("hyp4_control") is False or pr.get("hyp4_structural") is False:
                    fd.append(f"{key} {m}: Hyp4 violated")
        if sc == "physical" or sc == "algebraic":
            cls = "R" if sc == "physical" else "H"
            n_om = sum(1 for r in results if r["config"]["cls"] == cls for pr in r["pairs"].values()
                       if pr["category"] == "ordering-meaningful")
            per_m = {m: sum(1 for r in results if r["config"]["cls"] == cls
                            if r["pairs"].get(m, {}).get("category") == "ordering-meaningful") for m in PURE_ORDER}
            exp = EXPECTED_PAIRS[set_name][cls]
            if n_om < 0.9 * exp or any(v < 10 for v in per_m.values()):
                di.append(f"coverage: {n_om} ordering-meaningful class-{cls} pairs (expected {exp}); per mutant {per_m}")
        if inv:
            dec = "INVALID"
        elif fw:
            dec = "FAIL-W"
        elif fd:
            dec = "FAIL-D"
        elif di:
            dec = "DESIGN INSUFFICIENT"
        else:
            dec = "PASS" if set_name == "HO" else "PASS (exploratory; XP can never yield PASS)"
        out["scopes"][sc] = dict(decision=dec, fail_w=fw, fail_d=fd, design_insufficient=di)
    return out


# ============================================================================================ determinism, outputs
def determinism(cfgs, mods, grid=GRID):
    """§4.7 item 7: first 5 configurations, double outputs over the grid recomputed and compared bitwise."""
    wz, wzm = mods[0], mods[1]
    ok = True
    for cfg in cfgs[:5]:
        for s in grid:
            hs = []
            for _ in range(2):
                di = double_inputs(cfg, s)
                lead = wz.lead_from_delta(di["dlead"])
                GA = wz.wz_endpoint_blocks(di["Delta"], lead["lam_minus_1"])
                arrs = [GA[k] for k in KEYS.values()] + [wz.wz_fisher_lee(GA, lead["sin_theta"], lead["open"])["S"]]
                for m in wzm.MUTANTS:
                    if wzm.on_path(m, cfg["N"]):
                        r = wzm.wz_instrumented(di["Delta"], lead["lam_minus_1"], mutant=m)
                        arrs += [r[k] for k in KEYS.values()]
                hs.append(_hash_arrays(arrs))
            ok = ok and hs[0] == hs[1]
    return ok


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (Fraction,)):
        return f"{x.numerator}/{x.denominator}"
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, mpmath.mpf):
        return float(x)
    return x


POINT_COLS = ["set", "idx", "cls", "N_y", "N", "s", "mutant", "designated", "D", "Dstar_RW", "Dstar_L", "Dstar_tr_RW",
              "eta", "nu", "tau_abs_m", "nu_m", "snr", "Mag", "agree_ot_RW", "agree_do", "agree_do_ok", "agree_dt",
              "agree_dt_ok", "first_dev_ok", "item3", "tau_end", "judgeable", "cond_M_RW", "A_max_e_glob", "A_dS",
              "B_max_e_glob", "B_dS"] + [f"{q}_{b}" for b in BLOCKS for q in ("D", "Dtr_RW", "Dtr_L", "Dor_RW")]
PAIR_COLS = ["set", "idx", "cls", "N_y", "N", "mutant", "category", "q", "mu_est", "p_hat", "p_hat_star", "interval", "se",
             "n_points", "span", "window_adequate", "s_cross", "s_cross_secondary", "regime", "phi1_exact",
             "phi1_double", "hyp1", "hyp2", "hyp4", "series_check"] + [f"C{k}_{d}" for k in range(13) for d in (100, 200)] \
    + [f"z{k}" for k in range(13)]


def write_outputs(set_name, results, decision, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, f"results_{set_name}.json"), "w") as f:
        json.dump(_jsonable(dict(set=set_name, decision=decision, configurations=results)), f, sort_keys=True, indent=1)
    with open(os.path.join(out_dir, f"points_{set_name}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, POINT_COLS)
        w.writeheader()
        for r in results:
            c = r["config"]
            for p in r["points"]:
                if p.get("aborted"):
                    continue
                for m, mr in p["mutants"].items():
                    if m.startswith("_") or mr["status"] != "computed":
                        continue
                    row = dict(set=set_name, idx=c["idx"], cls=c["cls"], N_y=c["ny"], N=c["N"], s=p["s"], mutant=m,
                               designated=mr["designated"], D=mr["D"], Dstar_RW=mr["Dstar_RW"], Dstar_L=mr["Dstar_L"],
                               Dstar_tr_RW=mr["Dstar_tr_RW"], eta=mr["eta"], nu=mr["nu"], tau_abs_m=mr["tau_abs_m"],
                               nu_m=mr["nu_m"], snr=mr["snr"], Mag=mr["Mag"],
                               agree_ot_RW=mr["agree_oracle_transcription"]["R_W"],
                               agree_do=mr["agree_double_oracle"]["value"], agree_do_ok=mr["agree_double_oracle"]["ok"],
                               agree_dt=mr["agree_double_transcription"]["value"],
                               agree_dt_ok=mr["agree_double_transcription"]["ok"], first_dev_ok=mr["first_deviation"]["ok"],
                               item3=mr["item3"], tau_end=p["tau_end"], judgeable=p["judgeable"], cond_M_RW=p["cond"]["R_W"],
                               A_max_e_glob=max(v["e_glob"] for v in p["A"]["blocks"].values()), A_dS=p["A"]["dS_max"],
                               B_max_e_glob=p["B"]["max_e_glob"], B_dS=p["B"]["dS_max"])
                    for b in BLOCKS:
                        for q in ("D", "Dtr_RW", "Dtr_L", "Dor_RW"):
                            row[f"{q}_{b}"] = mr["blocks"][b][q]
                    w.writerow(row)
    with open(os.path.join(out_dir, f"pairs_{set_name}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, PAIR_COLS)
        w.writeheader()
        for r in results:
            c = r["config"]
            for m, pr in r["pairs"].items():
                reg = pr.get("regression", {})
                row = dict(set=set_name, idx=c["idx"], cls=c["cls"], N_y=c["ny"], N=c["N"], mutant=m,
                           category=pr["category"], q=pr.get("q"), mu_est=pr.get("mu_est"), p_hat=reg.get("p_hat"),
                           p_hat_star=reg.get("p_hat_star"), interval=reg.get("interval"), se=reg.get("se"),
                           n_points=reg.get("n_points"), span=reg.get("span"), window_adequate=pr.get("window_adequate"),
                           s_cross=pr.get("s_cross"), s_cross_secondary=pr.get("s_cross_secondary"),
                           regime=json.dumps(pr.get("regime")), phi1_exact=pr.get("phi1_exact"),
                           phi1_double=pr.get("phi1_double"), hyp1=pr.get("hyp1"), hyp2=pr.get("hyp2"),
                           hyp4=pr.get("hyp4_control", pr.get("hyp4_structural")),
                           series_check=(pr.get("series_check") or {}).get("status"))
                for k in range(13):
                    row[f"C{k}_100"] = (pr.get("Cnorm100") or [None] * 13)[k]
                    row[f"C{k}_200"] = (pr.get("Cnorm200") or [None] * 13)[k]
                    row[f"z{k}"] = (pr.get("tests") or [None] * 13)[k]
                w.writerow(row)


def write_report(decisions, out_dir):
    lines = ["# Coupling-scaling experiment v2 — report", ""]
    for d in decisions:
        lines.append(f"## Set {d['set']} ({'decides' if d['decides'] else 'exploratory; reported only'})")
        lines.append(f"- INVALID conditions: {len(d['invalid'])}")
        for sc, v in d["scopes"].items():
            title = {"physical": "class R + class D controls (physical scope)", "algebraic": "class H (algebraic stress test)"}[sc]
            lines.append(f"### {title}: **{v['decision']}**")
            for k in ("fail_w", "fail_d", "design_insufficient"):
                lines.append(f"- {k}: {len(v[k])}")
                lines += [f"  - {x}" for x in v[k][:50]]
        lines.append("")
    with open(os.path.join(out_dir, "coupling_scaling_v2_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="coupling-scaling v2 driver (gated)")
    ap.add_argument("--set", choices=("XP", "HO"), required=True)
    ap.add_argument("--out", default=OUT_DIR)
    a = ap.parse_args(argv)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    gate(a.out)                                              # raises GateRefused
    with open(os.path.join(a.out, f"manifest_{a.set}.json")) as f:
        man = json.load(f)
    cfgs = [config_from_manifest(r) for r in man["configurations"]]
    mods = _solvers()
    results = [run_configuration(c, GRID, mods) for c in cfgs]
    det = determinism(cfgs, mods)
    dec = decide(a.set, results, det)
    write_outputs(a.set, results, dec, a.out)
    write_report([dec], a.out)
    gate(a.out)                                              # §8.3: previous files and manifest unchanged after the run
    return 0


if __name__ == "__main__":
    sys.exit(main())
