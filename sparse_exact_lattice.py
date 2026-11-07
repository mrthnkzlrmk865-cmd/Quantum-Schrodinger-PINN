r"""
Tier 3 of B3Y: exact-lattice reference (docs/RGF_B3Y_VALIDATION.md, §I.6-I.8).

Computes the S-matrix of the mode-space lattice problem of spec §2.7 in EXTENDED precision
(np.longdouble, >= 64 mantissa bits) and keeps three objects apart:

  1. the exact lattice operator A  (exact inputs: the double E, E_perp and potential samples,
     dx = 2^-p, hence t = 1/(2 dx^2) exact; the DST-I definition of V_j and the lead
     elimination Sigma = -t Lambda are evaluated in extended precision);
  2. its floating-point representation fl(A), the merged double matrix a method assembles;
  3. the numerical solution of fl(A) x = b by a double-precision LU.

UNMERGED OPERATOR. With P_j = diag(E - E_perp) - V_j the exact operator acts as
    (A c)_j = P_j c_j + t (c_{j-1} - 2 c_j + c_{j+1}),   ghosts  c_0 = Lambda c_1,  c_{N+1} = Lambda c_N.
The large 2t is never added to the O(10^2) physical part. A double-precision SuperLU of the merged
matrix is used ONLY as a preconditioner: residuals are evaluated with the unmerged operator in extended
precision and the iterates are stored in extended precision (iterative refinement). A small residual of
the ROUNDED matrix is never used as evidence (it would measure the solver only).

Platform: np.longdouble must give >= 64 mantissa bits (x86-64 80-bit extended). Otherwise every entry
point raises ExtendedPrecisionUnavailable ("NOT EVALUABLE"); ordinary double precision is never substituted.

Independent of the RGF and of the other references: imports only the standard library, NumPy, SciPy and
mpmath (it does not import stationary, smatrix_spec, sparse_reference, sparse_wavefunction or any
validate_* module).

Contents
  platform_precision / require_extended_precision      gate X0
  lead_data_extended / lead_data_mpmath                 gate X1 (cancellation-free formulas vs mpmath)
  potential_blocks_extended                             cosine-sum DST-I blocks, extended precision
  UnmergedOperator, apply_matrix_extended               extended-precision operator applications
  iterative_refinement, refinement_gate                 gate X4 (contraction, iterations, operator residual)
  exact_lattice_S[_from_blocks]                         S from the scattering state AND from Fisher-Lee (gate X5),
                                                        preconditioner-ordering stability (gate X6), indicator e_ind (X7)
  solve_matrix_extended                                 exact solution of a ROUNDED double matrix (for the decomposition)
  error_decomposition, classify_*                       telescoping decomposition and the predeclared classification rules
  mp_single_channel_amplitudes, mp_dense_S              mpmath scalar references (X2, X3 and dense-mixing checks)
"""

from __future__ import annotations

import math
import platform
import sys
import time

import mpmath as mp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

LD = np.longdouble
CLD = np.clongdouble

# ----------------------------------------------------------------------
# predeclared values (docs/rgf_b3y_protocol.json; checked against it by the tests)
# ----------------------------------------------------------------------
MIN_MANTISSA_BITS = 64
LEAD_REL_TOL = 1e-17
SCALAR_ABS_TOL = 1e-13
CONTRACTION_MAX = 0.1
MAX_ITERATIONS = 10
ETA_MAX = 1e-17
FL_WF_TOL = 1e-13
ORDERING_TOL = 1e-13
ERROR_INDICATOR_MAX = 1e-11
VALIDATION_ABS = 1e-10
MPMATH_DPS = 50
GROWTH_PER_HALVING_MAX = 10.0
DOMINANCE_FACTOR = 2.0
TERM_RESOLUTION_FACTOR = 10.0
S_TARGET = 1e-10
THRESHOLD_TOL = 1e-12


class ExtendedPrecisionUnavailable(RuntimeError):
    """np.longdouble does not provide >= 64 mantissa bits: Tier 3 is NOT EVALUABLE on this platform."""


# ----------------------------------------------------------------------
# X0: platform precision
# ----------------------------------------------------------------------
def platform_precision() -> dict:
    fi = np.finfo(np.longdouble)
    bits = int(fi.nmant) + 1
    ok = bits >= MIN_MANTISSA_BITS
    return dict(longdouble_mantissa_bits=bits, longdouble_eps=float(fi.eps), longdouble_itemsize_bytes=int(np.dtype(np.longdouble).itemsize),
                float64_eps=float(np.finfo(np.float64).eps), machine=platform.machine(), system=platform.system(), sys_platform=sys.platform,
                numpy=np.__version__, python=platform.python_version(), required_bits=MIN_MANTISSA_BITS, evaluable=bool(ok),
                reason="ok" if ok else (f"np.longdouble has only {bits} mantissa bits (< {MIN_MANTISSA_BITS}): on {platform.system()}/{platform.machine()} "
                                        f"longdouble is not an extended-precision type"))


def require_extended_precision() -> dict:
    info = platform_precision()
    if not info["evaluable"]:
        raise ExtendedPrecisionUnavailable("NOT EVALUABLE: " + info["reason"] + ". Tier 3 stops; ordinary double precision is never substituted.")
    return info


def eps_ld() -> float:
    return float(np.finfo(np.longdouble).eps)


def hopping(dx: float) -> LD:
    d = LD(dx)
    return LD(1) / (LD(2) * d * d)


# ----------------------------------------------------------------------
# X1: lead data (cancellation-free, extended) and the mpmath scalar reference
# ----------------------------------------------------------------------
def lead_data_extended(E: float, E_perp, dx: float) -> dict:
    """x_n = (E - E_perp_n) dx^2 is formed directly (no c = 1 - x).
    Open (0 < x < 2): lambda = (1 - x) + i sqrt(x(2 - x)), Gamma = 2 t sqrt(x(2 - x)).
    Closed x < 0: lambda = 1/(c + sqrt(y(2 + y))), y = -x, c = 1 + y.  x > 2 (lattice artifact): lambda = -1/(z + sqrt(x(x-2))), z = x - 1."""
    require_extended_precision()
    Ep = np.asarray(E_perp, dtype=np.float64).astype(LD)
    d = LD(dx)
    x = (LD(E) - Ep) * d * d
    if np.any(np.abs(x) < THRESHOLD_TOL) or np.any(np.abs(x - 2) < THRESHOLD_TOL):
        raise ValueError("energy at a channel threshold or band edge")
    op = (x > 0) & (x < 2)
    lam = np.empty(x.shape, dtype=CLD)
    gam = np.zeros(x.shape, dtype=LD)
    t = hopping(dx)
    xo = x[op]
    s = np.sqrt(xo * (2 - xo))
    lam[op] = (1 - xo).astype(CLD) + s.astype(CLD) * CLD(1j)
    gam[op] = 2 * t * s
    below = x < 0
    y = -x[below]
    lam[below] = (1 / ((1 + y) + np.sqrt(y * (2 + y)))).astype(CLD)
    above = x > 2
    xa = x[above]
    lam[above] = (-1 / ((xa - 1) + np.sqrt(xa * (xa - 2)))).astype(CLD)
    return dict(x=x, open=op, lam=lam, gamma=gam, t=t, sin_theta=np.where(op, np.sqrt(np.abs(x * (2 - x))), 0))


def lead_data_mpmath(E: float, E_perp, dx: float, dps: int = MPMATH_DPS) -> dict:
    """Scalar reference with different formulas: c = 1 - x exactly (50 digits), theta = acos(c), lambda = e^{i theta}, Gamma = 2 t sin(theta);
    closed: the root of lambda^2 - 2 c lambda + 1 with |lambda| < 1."""
    lam, gam, op = [], [], []
    with mp.workdps(dps):
        d = mp.mpf(float(dx))
        t = 1 / (2 * d * d)
        for ep in np.asarray(E_perp, dtype=np.float64):
            x = (mp.mpf(float(E)) - mp.mpf(float(ep))) * d * d
            c = 1 - x
            if abs(c) < 1:
                th = mp.acos(c)
                lam.append(mp.expj(th)); gam.append(2 * t * mp.sin(th)); op.append(True)
            else:
                r = mp.sqrt(c * c - 1)
                l1, l2 = c - r, c + r
                lam.append(l1 if abs(l1) < abs(l2) else l2); gam.append(mp.mpf(0)); op.append(False)
    return dict(lam=lam, gamma=gam, open=np.array(op))


def _to_mp(z, dps: int = MPMATH_DPS):
    f = lambda v: mp.mpf(np.format_float_positional(np.longdouble(v), precision=24, unique=False))
    z = np.clongdouble(z)
    return mp.mpc(f(z.real), f(z.imag))


def lead_data_vs_mpmath(E: float, E_perp, dx: float) -> dict:
    """Gate X1: relative differences of lambda and Gamma (every channel) against mpmath."""
    ext, ref = lead_data_extended(E, E_perp, dx), lead_data_mpmath(E, E_perp, dx)
    assert np.array_equal(ext["open"], ref["open"])
    rel_lam, rel_gam = 0.0, 0.0
    with mp.workdps(MPMATH_DPS):
        for i in range(len(ref["lam"])):
            rel_lam = max(rel_lam, float(abs(_to_mp(ext["lam"][i]) - ref["lam"][i]) / abs(ref["lam"][i])))
            if ref["open"][i]:
                g = mp.mpf(np.format_float_positional(ext["gamma"][i], precision=24, unique=False))
                rel_gam = max(rel_gam, float(abs(g - ref["gamma"][i]) / ref["gamma"][i]))
    return dict(rel_lambda=rel_lam, rel_gamma=rel_gam, passed=bool(max(rel_lam, rel_gam) <= LEAD_REL_TOL))


# ----------------------------------------------------------------------
# potential blocks (extended precision; cosine sums of the DST-I quadrature)
# ----------------------------------------------------------------------
def potential_blocks_extended(V_samples, n_modes: int, n_quad: int) -> np.ndarray:
    """V_nm = (g_{|n-m|} - g_{n+m})/(N_q+1), g_k = sum_l V(x_j, y_l) cos(pi k l/(N_q+1)); the cosine argument is reduced with exact
    integer arithmetic (k l mod 2(N_q+1)) so no accuracy is lost to the size of the argument."""
    require_extended_precision()
    if n_quad < 2 * n_modes:
        raise ValueError("n_quad must be >= 2 n_modes")
    V = np.asarray(V_samples, dtype=np.float64).astype(LD)
    if V.ndim != 2 or V.shape[1] != n_quad:
        raise ValueError(f"V_samples must have shape (N, {n_quad})")
    m = np.outer(np.arange(0, 2 * n_modes + 1), np.arange(1, n_quad + 1)) % (2 * (n_quad + 1))
    pi = np.arccos(LD(-1))
    cosm = np.cos(pi * m.astype(LD) / LD(n_quad + 1))
    g = V @ cosm.T
    n = np.arange(1, n_modes + 1)
    return (g[:, np.abs(n[:, None] - n[None, :])] - g[:, n[:, None] + n[None, :]]) / LD(n_quad + 1)


# ----------------------------------------------------------------------
# operators in extended precision
# ----------------------------------------------------------------------
class UnmergedOperator:
    """A c = P c + t (c_{j-1} - 2 c_j + c_{j+1}) with lead ghosts; P_j = diag(E - E_perp) - V_j kept separate from 2t."""

    def __init__(self, E, E_perp, V_blocks, dx, lam):
        require_extended_precision()
        Vb = np.asarray(V_blocks).astype(LD)
        self.N, self.ny = Vb.shape[0], Vb.shape[1]
        Ep = np.asarray(E_perp, dtype=np.float64).astype(LD)
        self.P = np.zeros_like(Vb)
        self.P -= Vb
        idx = np.arange(self.ny)
        self.P[:, idx, idx] += LD(E) - Ep
        self.t = hopping(dx)
        self.lam = np.asarray(lam, dtype=CLD)

    def apply(self, X: np.ndarray) -> np.ndarray:
        Y = np.empty(X.shape, dtype=CLD)
        Y.real = np.matmul(self.P, np.ascontiguousarray(X.real))
        Y.imag = np.matmul(self.P, np.ascontiguousarray(X.imag))
        lam = self.lam[None, :, None]
        left = np.empty_like(X); right = np.empty_like(X)
        left[1:] = X[:-1]; left[0] = lam[0] * X[0]
        right[:-1] = X[1:]; right[-1] = lam[0] * X[-1]
        Y += self.t * (left - 2 * X + right)
        return Y


def apply_matrix_extended(M: sp.spmatrix, X: np.ndarray, chunk_rows: int = 20000) -> np.ndarray:
    """y = M x for a ROUNDED double-precision sparse matrix, evaluated exactly in extended precision (entries are exact doubles).
    X: (n, nrhs) clongdouble."""
    M = M.tocsr()
    n, nr = X.shape
    Y = np.zeros((n, nr), dtype=CLD)
    indptr, indices = M.indptr, M.indices
    data = M.data.astype(np.complex128).astype(CLD)
    for r0 in range(0, n, chunk_rows):
        r1 = min(r0 + chunk_rows, n)
        a, b = indptr[r0], indptr[r1]
        if b == a:
            continue
        prod = data[a:b, None] * X[indices[a:b]]
        starts = indptr[r0:r1] - a
        nonempty = indptr[r0 + 1:r1 + 1] > indptr[r0:r1]
        sums = np.add.reduceat(prod, starts[nonempty], axis=0)
        Y[r0:r1][nonempty] = sums
    return Y


# ----------------------------------------------------------------------
# merged double matrix: PRECONDITIONER ONLY
# ----------------------------------------------------------------------
def assemble_preconditioner(E, E_perp, V_blocks, dx, lam) -> sp.csc_matrix:
    """Merged double-precision matrix, (E - 2t) - E_perp - V on the diagonal blocks, +t Lambda at the two ends, +t off the diagonal."""
    Vd = np.asarray(V_blocks).astype(np.float64)
    N, ny = Vd.shape[:2]
    t = float(hopping(dx))
    lam_d = np.asarray(lam).astype(np.complex128)
    base = (float(E) - 2.0 * t) * np.eye(ny) - np.diag(np.asarray(E_perp, dtype=np.float64))
    blocks = [(base - Vd[j]).astype(complex) for j in range(N)]
    blocks[0] = blocks[0] + t * np.diag(lam_d)
    blocks[-1] = blocks[-1] + t * np.diag(lam_d)
    D = sp.block_diag([sp.csr_matrix(b) for b in blocks], format="csr")
    hop = sp.diags([np.ones(N - 1), np.ones(N - 1)], [-1, 1], shape=(N, N))
    return (D + t * sp.kron(hop, sp.identity(ny), format="csr")).tocsc()


def condition_estimate_1norm(M: sp.csc_matrix, lu) -> float:
    n = M.shape[0]
    lo = spla.LinearOperator((n, n), matvec=lambda v: lu.solve(v), rmatvec=lambda v: lu.solve(v, trans="H"),
                             matmat=lambda m: lu.solve(m), rmatmat=lambda m: lu.solve(m, trans="H"), dtype=complex)
    return float(spla.norm(M, 1) * spla.onenormest(lo, t=4))


# ----------------------------------------------------------------------
# sources and S extraction (extended precision)
# ----------------------------------------------------------------------
def build_sources(lead: dict, N: int, ny: int) -> np.ndarray:
    """Unit incident flux a = Gamma^{-1/2}; source s = -t (1/lambda - lambda) a at slice 1 (left) or N (right)."""
    o = np.where(lead["open"])[0]
    no = o.size
    a = lead["gamma"][o] ** LD(-0.5)
    src = (-lead["t"] * (1 / lead["lam"][o] - lead["lam"][o]) * a.astype(CLD)).astype(CLD)
    B = np.zeros((N, ny, 2 * no), dtype=CLD)
    B[0, o, np.arange(no)] = src
    B[N - 1, o, no + np.arange(no)] = src
    return B


def S_from_state(X: np.ndarray, lead: dict) -> dict:
    o = np.where(lead["open"])[0]
    no = o.size
    sg = np.sqrt(lead["gamma"][o]).astype(CLD)[:, None]
    a = (lead["gamma"][o] ** LD(-0.5)).astype(CLD)
    c1, cN = X[0], X[-1]
    r_L = sg * (c1[o][:, :no] - np.diag(a)); t_RL = sg * cN[o][:, :no]
    t_LR = sg * c1[o][:, no:]; r_R = sg * (cN[o][:, no:] - np.diag(a))
    return dict(S=np.block([[r_L, t_LR], [t_RL, r_R]]), r_L=r_L, t_RL=t_RL, t_LR=t_LR, r_R=r_R, n_open=int(no))


def S_from_green_columns(X: np.ndarray, lead: dict) -> dict:
    """Fisher-Lee on Green columns X = G e (unit right-hand sides): t_RL = i G^{1/2}_R G_N1 G^{1/2}_L, r_L = i G^{1/2} G_11 G^{1/2} - I, ..."""
    o = np.where(lead["open"])[0]
    no = o.size
    sg = np.sqrt(lead["gamma"][o]).astype(CLD)
    I = np.eye(no, dtype=CLD)
    G11, GN1, G1N, GNN = X[0][o][:, :no], X[-1][o][:, :no], X[0][o][:, no:], X[-1][o][:, no:]
    f = lambda G: CLD(1j) * sg[:, None] * G * sg[None, :]
    r_L, t_RL, t_LR, r_R = f(G11) - I, f(GN1), f(G1N), f(GNN) - I
    return dict(S=np.block([[r_L, t_LR], [t_RL, r_R]]), r_L=r_L, t_RL=t_RL, t_LR=t_LR, r_R=r_R, n_open=int(no))


def build_unit_sources(lead: dict, N: int, ny: int) -> np.ndarray:
    o = np.where(lead["open"])[0]
    no = o.size
    B = np.zeros((N, ny, 2 * no), dtype=CLD)
    B[0, o, np.arange(no)] = 1
    B[N - 1, o, no + np.arange(no)] = 1
    return B


# ----------------------------------------------------------------------
# iterative refinement (extended-precision operator residual) and gate X4
# ----------------------------------------------------------------------
def iterative_refinement(apply_op, solve_prec, B: np.ndarray, norm_A1: float, extract=None, max_iter: int = MAX_ITERATIONS) -> dict:
    """x_0 = P^{-1} b; r_k = b - A x_k (A: unmerged operator, extended precision); d_k = P^{-1} r_k (double); x_{k+1} = x_k + d_k.
    P^{-1} is the double LU. Iterates are stored in extended precision. `extract(x)` -> S (for the per-iterate S difference)."""
    e = eps_ld()
    X = solve_prec(B)
    hist, S_prev, d_prev = [], None, None
    bn = float(np.max(np.abs(B)))
    converged = False
    for k in range(max_iter):
        R = B - apply_op(X)
        rn = float(np.max(np.abs(R)))
        D = solve_prec(R)
        dn = float(np.max(np.abs(D)))
        X = X + D
        xn = float(np.max(np.abs(X)))
        rec = dict(k=k, residual_inf=rn, correction_inf=dn, x_inf=xn, eta=rn / (bn + norm_A1 * xn),
                   contraction=(dn / d_prev) if d_prev else float("nan"))
        if extract is not None:
            S = extract(X)
            rec["dS_vs_previous"] = float(np.max(np.abs(S - S_prev))) if S_prev is not None else float("nan")
            S_prev = S
        hist.append(rec)
        d_prev = dn
        if k >= 1 and rec["eta"] <= ETA_MAX:                    # predeclared operator-residual target reached (>= 2 iterations, so an iterate difference exists)
            converged = True
            break
    Rf = B - apply_op(X)
    rn = float(np.max(np.abs(Rf)))
    return dict(X=X, history=hist, converged=converged, residual_inf_final=rn,
                eta_final=rn / (bn + norm_A1 * float(np.max(np.abs(X)))), n_iterations=len(hist))


def refinement_gate(result: dict) -> dict:
    """Gate X4: converged (relative operator residual eta <= 1e-17 reached), <= 10 iterations, and contraction <= 0.1 for every step taken while the
    previous relative operator residual still exceeded the 1e-17 target (below it the iteration is at the arithmetic floor, where the correction
    stagnates by construction)."""
    h = result["history"]
    worst = 0.0
    for k in range(1, len(h)):
        if h[k - 1]["eta"] > ETA_MAX:
            worst = max(worst, h[k]["contraction"])
    ok = bool(result["converged"] and result["n_iterations"] <= MAX_ITERATIONS and worst <= CONTRACTION_MAX and result["eta_final"] <= ETA_MAX)
    return dict(passed=ok, worst_contraction=worst, n_iterations=result["n_iterations"], eta_final=result["eta_final"], converged=result["converged"])


# ----------------------------------------------------------------------
# the exact-lattice S-matrix
# ----------------------------------------------------------------------
def exact_lattice_S_from_blocks(E, E_perp, V_blocks, dx, *, permc_spec: str = "MMD_AT_PLUS_A", max_iter: int = MAX_ITERATIONS,
                                with_fisher_lee: bool = True, independent_green_solve: bool = True, condition: bool = True,
                                return_solution: bool = False) -> dict:
    prec = require_extended_precision()
    t0 = time.time()
    lead = lead_data_extended(E, E_perp, dx)
    Vb = np.asarray(V_blocks)
    N, ny = Vb.shape[0], Vb.shape[1]
    op = UnmergedOperator(E, E_perp, Vb, dx, lead["lam"])
    Mp = assemble_preconditioner(E, E_perp, Vb, dx, lead["lam"])
    lu = spla.splu(Mp, permc_spec=permc_spec)
    norm1 = float(spla.norm(Mp, 1))
    kappa = condition_estimate_1norm(Mp, lu) if condition else float("nan")
    solve_prec = lambda R: lu.solve(R.astype(np.complex128).reshape(N * ny, -1)).reshape(N, ny, -1).astype(CLD)
    B = build_sources(lead, N, ny)
    ref = iterative_refinement(op.apply, solve_prec, B, norm1, extract=lambda X: S_from_state(X, lead)["S"], max_iter=max_iter)
    st = S_from_state(ref["X"], lead)
    hist = ref["history"]
    d_iter = hist[-1].get("dS_vs_previous", float("nan")) if len(hist) > 1 else float("nan")
    e_ind = max(d_iter if np.isfinite(d_iter) else 0.0, kappa * ref["eta_final"] if np.isfinite(kappa) else 0.0)
    out = dict(S=st["S"], r_L=st["r_L"], t_RL=st["t_RL"], t_LR=st["t_LR"], r_R=st["r_R"], n_open=st["n_open"], lead=lead,
               history=hist, converged=ref["converged"], eta_final=ref["eta_final"], residual_inf_final=ref["residual_inf_final"],
               refinement_gate=refinement_gate(ref), kappa1=kappa, error_indicator=e_ind,
               error_indicator_note="convergence/error indicator, not a rigorous uncertainty bound",
               permc_spec=permc_spec, preconditioner_lu_nnz=int(lu.L.nnz + lu.U.nnz), platform=prec, n_slices=N, n_modes=ny)
    if with_fisher_lee:
        if independent_green_solve:
            Bu = build_unit_sources(lead, N, ny)
            rg = iterative_refinement(op.apply, solve_prec, Bu, norm1, max_iter=max_iter)
            fl = S_from_green_columns(rg["X"], lead)
            out["green_refinement_gate"] = refinement_gate(rg)
        else:
            Xu = ref["X"] / B.sum(axis=(0, 1))[None, None, :]            # each column of B has exactly one non-zero entry (the source)
            fl = S_from_green_columns(Xu, lead)
        out["S_fisher_lee"] = fl["S"]
        out["fl_vs_state_max_abs"] = float(np.max(np.abs(fl["S"] - st["S"])))
    if return_solution:
        out["X"] = ref["X"]
    out["wall_time_s"] = time.time() - t0
    return out


def exact_lattice_S(E, E_perp, V_samples, dx, n_modes: int, n_quad: int, **kw) -> dict:
    """V_samples: (N, N_q) double samples of V(x_j, y_l); the extended-precision blocks are formed here (cosine sums)."""
    Vb = potential_blocks_extended(V_samples, n_modes, n_quad)
    out = exact_lattice_S_from_blocks(E, E_perp, Vb, dx, **kw)
    out["V_blocks"] = Vb
    return out


def ordering_stability(E, E_perp, V_blocks, dx, orderings=("COLAMD", "MMD_AT_PLUS_A")) -> dict:
    """Gate X6: the extended-precision S must not depend on the preconditioner ordering."""
    res = [exact_lattice_S_from_blocks(E, E_perp, V_blocks, dx, permc_spec=o, with_fisher_lee=False, condition=False) for o in orderings]
    d = float(max(np.max(np.abs(res[0]["S"] - r["S"])) for r in res[1:]))
    return dict(orderings=list(orderings), max_abs_difference=d, passed=bool(d <= ORDERING_TOL))


# ----------------------------------------------------------------------
# exact solution of a ROUNDED matrix (for the representation / solver decomposition)
# ----------------------------------------------------------------------
def solve_matrix_extended(M: sp.spmatrix, B: np.ndarray, *, permc_spec: str = "MMD_AT_PLUS_A", max_iter: int = MAX_ITERATIONS,
                          corner_swap: dict | None = None) -> dict:
    """Solve M x = b to extended precision, where M is a rounded double sparse matrix (its entries are exact doubles). B: (n, nrhs) clongdouble.
    The residual is evaluated EXACTLY (extended precision) from the stored entries; the double LU is the preconditioner. If corner_swap is given
    ({'t', 'ny', 'N', 'lam_double', 'lam_ext'}) the double lead terms t*Lambda_d in the two corner blocks are replaced by the exact extended
    t*Lambda: the 'hybrid' operator (rounded bulk, exact lead) used to isolate bulk-operator rounding."""
    require_extended_precision()
    Mc = M.tocsc()
    lu = spla.splu(Mc, permc_spec=permc_spec)
    n = Mc.shape[0]
    norm1 = float(spla.norm(Mc, 1))
    solve_prec = lambda R: lu.solve(R.astype(np.complex128)).astype(CLD)

    def apply(X):
        Y = apply_matrix_extended(Mc, X)
        if corner_swap is not None:
            ny, N = corner_swap["ny"], corner_swap["N"]
            dd = (corner_swap["t"] * np.asarray(corner_swap["lam_double"], dtype=np.complex128)).astype(CLD)
            de = (hopping_from_t(corner_swap["t"]) * np.asarray(corner_swap["lam_ext"], dtype=CLD)).astype(CLD)
            for sl in (slice(0, ny), slice((N - 1) * ny, N * ny)):
                Y[sl] += (de - dd)[:, None] * X[sl]
        return Y

    ref = iterative_refinement(apply, solve_prec, B.astype(CLD), norm1, max_iter=max_iter)
    ref["refinement_gate"] = refinement_gate(ref)
    return ref


def hopping_from_t(t) -> LD:
    return LD(t)


def error_decomposition(S_double, S_exact, S_fl, S_hybrid=None) -> dict:
    """Telescoping identity  S_double - S_exact = (S_hybrid - S_exact) + (S_fl - S_hybrid) + (S_double - S_fl):
    bulk-operator rounding + lead data + linear-solver/extraction rounding.  Without S_hybrid: representation (S_fl - S_exact) + solver."""
    n = lambda A: float(np.max(np.abs(A)))
    out = dict(total=n(S_double - S_exact), representation_total=n(S_fl - S_exact), solver=n(S_double - S_fl))
    if S_hybrid is not None:
        out.update(bulk_rounding=n(S_hybrid - S_exact), lead_data=n(S_fl - S_hybrid))
    out["identity_residual"] = n((S_double - S_exact) - ((S_fl - S_exact) + (S_double - S_fl)))
    return out


# ----------------------------------------------------------------------
# predeclared classification rules (docs/RGF_B3Y_VALIDATION.md §I.8.3, revision 1)
# ----------------------------------------------------------------------
def precision_consistency(eps_by_dx: dict) -> dict:
    """OPERATIONAL classification (not a proof of round-off dominance): eps(1/512) <= 1e-10 and growth <= 10 per halving of dx.
    eps_by_dx maps inverse dx (512, 1024, 2048) to eps_m; missing levels are skipped."""
    levels = sorted(eps_by_dx)
    ok0 = bool(512 in eps_by_dx and eps_by_dx[512] <= S_TARGET)
    growth = [eps_by_dx[b] / eps_by_dx[a] for a, b in zip(levels[:-1], levels[1:]) if eps_by_dx[a] > 0]
    return dict(eps_512_within_target=ok0, growth_per_halving=growth, passed=bool(ok0 and all(g <= GROWTH_PER_HALVING_MAX for g in growth)))


def representation_dominance(terms: dict | None, e_ind: float) -> dict:
    """Supported only if (M2/M3, decomposition available) bulk_rounding >= 2 x max(lead_data, solver) and every term > 10 x e_ind."""
    if terms is None or any(k not in terms for k in ("bulk_rounding", "lead_data", "solver")):
        return dict(available=False, supported=False, resolved=False, reason="no explicit matrix: only the total error is available")
    resolved = all(terms[k] > TERM_RESOLUTION_FACTOR * e_ind for k in ("bulk_rounding", "lead_data", "solver"))
    dominant = terms["bulk_rounding"] >= DOMINANCE_FACTOR * max(terms["lead_data"], terms["solver"])
    return dict(available=True, resolved=bool(resolved), dominant=bool(dominant), supported=bool(resolved and dominant),
                reason="supported" if (resolved and dominant) else ("terms unresolved" if not resolved else "bulk-rounding term does not dominate"))


def classify_exceedance(eps_by_dx: dict, terms: dict | None, e_ind: float, claim_representation_dominance: bool) -> str:
    """Classification of a method whose error against the exact lattice exceeds 1e-10 (never PASS)."""
    pc = precision_consistency(eps_by_dx)
    if not pc["passed"]:
        return "UNEXPLAINED"
    if not claim_representation_dominance:
        return "PRECISION-CONSISTENT"
    return "PRECISION-LIMITED (measured; representation-dominated)" if representation_dominance(terms, e_ind)["supported"] else "UNEXPLAINED"


# ----------------------------------------------------------------------
# mpmath scalar references
# ----------------------------------------------------------------------
def mp_single_channel_amplitudes(E: float, E_perp_n: float, V_x, dx: float, dps: int = MPMATH_DPS):
    """(r, t) of one decoupled channel on the lattice, left incidence, convention P: exact backward transfer recursion from the exact inputs."""
    with mp.workdps(dps):
        d = mp.mpf(float(dx)); t = 1 / (2 * d * d)
        x = (mp.mpf(float(E)) - mp.mpf(float(E_perp_n))) * d * d
        th = mp.acos(1 - x)
        lam = mp.expj(th)
        e0, ep = mp.mpf(float(E)), mp.mpf(float(E_perp_n))
        nxt, cur = lam, mp.mpf(1)
        for j in range(len(V_x), 0, -1):
            prev = ((2 * t + ep + mp.mpf(float(V_x[j - 1])) - e0) * cur - t * nxt) / t
            nxt, cur = cur, prev
        psi1, psi0 = nxt, cur
        A = mp.matrix([[1, 1], [1 / lam, lam]])
        sol = mp.lu_solve(A, mp.matrix([psi1, psi0]))
        a, b = sol[0], sol[1]
        return b / a, 1 / a


def mp_dense_S(E: float, E_perp, V_blocks, dx: float, dps: int = MPMATH_DPS) -> np.ndarray:
    """Exact S (mpmath dense solve of the unmerged operator with the exact leads) for a SMALL device with dense (mixing) blocks; returns an object array."""
    with mp.workdps(dps):
        Vb = np.asarray(V_blocks, dtype=np.float64)
        N, ny = Vb.shape[:2]
        d = mp.mpf(float(dx)); t = 1 / (2 * d * d)
        lam, gam, op = [], [], []
        for ep in np.asarray(E_perp, dtype=np.float64):
            x = (mp.mpf(float(E)) - mp.mpf(float(ep))) * d * d
            c = 1 - x
            if abs(c) < 1:
                th = mp.acos(c); lam.append(mp.expj(th)); gam.append(2 * t * mp.sin(th)); op.append(True)
            else:
                r = mp.sqrt(c * c - 1); l1, l2 = c - r, c + r
                lam.append(l1 if abs(l1) < abs(l2) else l2); gam.append(mp.mpf(0)); op.append(False)
        o = [i for i in range(ny) if op[i]]
        no = len(o)
        n = N * ny
        A = mp.zeros(n, n)
        for j in range(N):
            for a in range(ny):
                for b in range(ny):
                    v = -mp.mpf(float(Vb[j, a, b]))
                    if a == b:
                        v += mp.mpf(float(E)) - mp.mpf(float(E_perp[a])) - 2 * t
                    A[j * ny + a, j * ny + b] += v
                if j == 0 or j == N - 1:
                    A[j * ny + a, j * ny + a] += t * lam[a]
                if j + 1 < N:
                    A[j * ny + a, (j + 1) * ny + a] += t
                if j > 0:
                    A[j * ny + a, (j - 1) * ny + a] += t
        B = mp.zeros(n, 2 * no)
        for i, ch in enumerate(o):
            a_i = gam[ch] ** mp.mpf(-0.5)
            s = -t * (1 / lam[ch] - lam[ch]) * a_i
            B[ch, i] = s
            B[(N - 1) * ny + ch, no + i] = s
        X = mp.zeros(n, 2 * no)
        for c in range(2 * no):                                  # mpmath's lu_solve takes one right-hand side (the LU is cached on A)
            X[:, c] = mp.lu_solve(A, B[:, c])
        S = np.empty((2 * no, 2 * no), dtype=object)
        for i in range(2 * no):
            for k in range(2 * no):
                sg = mp.sqrt(gam[o[i % no]])
                if i < no:
                    val = sg * (X[o[i], k] - (gam[o[k]] ** mp.mpf(-0.5) if k == i else 0))      # c_1 rows (left outgoing)
                    if k >= no:
                        val = sg * X[o[i], k]
                else:
                    val = sg * X[(N - 1) * ny + o[i - no], k]
                    if k >= no and (k - no) == (i - no):
                        val = sg * (X[(N - 1) * ny + o[i - no], k] - gam[o[k - no]] ** mp.mpf(-0.5))
                S[i, k] = val
        return S
