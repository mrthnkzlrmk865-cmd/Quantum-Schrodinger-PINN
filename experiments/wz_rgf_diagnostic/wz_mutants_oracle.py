r"""Specification-derived mutant oracle of §4.2 (PREDICTIONS_coupling_scaling_v2.md, Revision 3.7).

For every mutant, the exact mutant value of its designated block X_m is computed from the direct prefix and full inverses
(Ĝ units) by the closed forms of §4.2. No recursion code is used (§4.1). The formulas are written against a small matrix
backend so that the same code runs on mpmath references (driver; `MpBackend`) and on exact rationals (tests).

Notation (§4.2): Ĝ^(j) = (M^(j)/t)^-1 (prefix, left lead only), Ĝ = (M/t)^-1 (full); P_j = -Ĝ^(j)_jj, P_1 = -Ĝ^(1)_11;
chi_j = Ĝ^(j)_j1, psi_j = Ĝ^(j)_1j, theta_j = psi_{j-1} P_j chi_{j-1}; Z^-1 = Ĝ_NN; sums over loop slices j = 2..N-1.
Blocks are 1-indexed slices of N_y x N_y.

Status of every mutant on a given N:
  * off_path      — N < MIN_N[m], or m = M11 (an N = 1 boundary mutant; the experiment has N >= 2). No value is returned:
                    an off-path mutant is never a zero discrepancy.
  * computed      — the designated block (and, for M14, its two secondary blocks named in §4.2) is returned.
  * Identical-by-algebra (M5 at N = 3, M8 at N = 2) and equivalent (M12) mutants are *computed* by their formulas; their
    category is assigned by the driver (§5.1), not here.
Only the designated block X_m is specified by §4.2 (plus M14's Ĝ_N1 -> -Ĝ_N1 and Ĝ_1N -> -Ĝ_1N); other blocks of the
mutant output are not specified by the oracle and are not returned.
"""
from __future__ import annotations

MUTANTS = ("M1", "M2", "M3", "M4", "M4L", "M4R", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14",
           "M15", "M16")
MIN_N = {"M1": 3, "M2": 3, "M3": 3, "M4": 3, "M4L": 3, "M4R": 3, "M5": 3, "M9": 3, "M13": 3,
         "M6": 2, "M7": 2, "M8": 2, "M10": 2, "M12": 2, "M14": 2, "M15": 2, "M16": 2}
DESIGNATED = {"M1": "G_N1", "M2": "G_1N", "M3": "G_N1", "M4": "G_11", "M4L": "G_11", "M4R": "G_11", "M5": "G_11",
              "M6": "G_N1", "M7": "G_1N", "M8": "G_11", "M9": "G_N1", "M10": "G_NN", "M12": "G_11", "M13": "G_11",
              "M14": "G_11", "M15": "G_N1", "M16": "G_11"}


class OffPath(ValueError):
    """The mutant is not on the code path for this N (no oracle value exists)."""


def on_path(m, N):
    if m not in MUTANTS:
        raise ValueError(f"unknown mutant {m!r}")
    if m == "M11":
        return N == 1
    return N >= MIN_N[m]


class MpBackend:
    """mpmath matrices in a given context (the driver's reference precision)."""

    def __init__(self, ctx):
        self.ctx = ctx

    def mul(self, A, B):
        return A * B

    def add(self, A, B):
        return A + B

    def sub(self, A, B):
        return A - B

    def neg(self, A):
        return -A

    def scal(self, c, A):
        return self.ctx.mpf(c) * A

    def inv(self, A):
        return self.ctx.inverse(A)

    def block(self, X, a, b, ny):
        return X[(a - 1) * ny:a * ny, (b - 1) * ny:b * ny]

    def diag(self, vals):
        return self.ctx.diag([self.ctx.mpc(v) for v in vals])


class _Blocks:
    def __init__(self, be, G_full, G_prefix, N, ny):
        self.be, self.Gf, self.Gp, self.N, self.ny = be, G_full, G_prefix, N, ny

    def full(self, a, b):
        return self.be.block(self.Gf, a, b, self.ny)

    def pre(self, j, a, b):
        if j not in self.Gp:
            raise KeyError(f"prefix inverse Ĝ^({j}) required but not supplied")
        return self.be.block(self.Gp[j], a, b, self.ny)

    def P(self, j):
        return self.be.neg(self.pre(j, j, j))            # P_j = -Ĝ^(j)_jj (P_1 = -Ĝ^(1)_11)

    def chi(self, j):
        return self.pre(j, j, 1)

    def psi(self, j):
        return self.pre(j, 1, j)

    def theta(self, j):
        be = self.be
        return be.mul(be.mul(self.psi(j - 1), self.P(j)), self.chi(j - 1))


def _prod(be, mats):
    X = mats[0]
    for Y in mats[1:]:
        X = be.mul(X, Y)
    return X


def oracle(m, N, ny, G_full, G_prefix, lam_minus_1=None, backend=None):
    """Exact mutant value of the designated block of mutant m (§4.2).

    G_full: (M/t)^-1 (n_tot x n_tot); G_prefix: {j: (M^(j)/t)^-1} for j = 1..N-1; lam_minus_1: diagonal of Lambda - I
    (needed by M10 only). Returns {"status": "computed", "designated": name, "X": block, "secondary": {...}}.
    Raises OffPath for an off-path mutant and ValueError for an unknown label."""
    be = backend
    if not on_path(m, N) or m == "M11":
        raise OffPath(f"{m} is off path at N = {N}")
    B = _Blocks(be, G_full, G_prefix, N, ny)
    GNN, GN1, G1N, G11 = B.full(N, N), B.full(N, 1), B.full(1, N), B.full(1, 1)
    loop = range(2, N)                                    # j = 2..N-1
    sec = {}
    if m in ("M1", "M3"):
        X = be.neg(_prod(be, [GNN, B.pre(1, 1, 1)] + [B.P(j) for j in loop]))
    elif m == "M2":
        X = be.neg(_prod(be, [B.P(j) for j in reversed(loop)] + [B.pre(1, 1, 1), GNN]))
    elif m in ("M4", "M4L", "M4R"):
        X = G11
        for j in loop:
            Pj = B.P(j)
            Pk = _prod(be, [Pj, Pj, Pj]) if m == "M4" else _prod(be, [Pj, Pj])
            X = be.sub(X, _prod(be, [B.psi(j - 1), be.sub(Pk, Pj), B.chi(j - 1)]))
    elif m == "M5":
        X = G11
        for j in loop:
            X = be.sub(X, be.sub(_prod(be, [B.chi(j - 1), B.P(j), B.psi(j - 1)]), B.theta(j)))
    elif m in ("M6", "M15"):
        X = be.neg(be.mul(B.chi(N - 1), GNN))
    elif m == "M7":
        X = be.neg(be.mul(GNN, B.psi(N - 1)))
    elif m == "M8":
        X = be.add(G11, be.sub(_prod(be, [B.chi(N - 1), GNN, B.psi(N - 1)]), _prod(be, [B.psi(N - 1), GNN, B.chi(N - 1)])))
    elif m == "M9":
        # -Z^-1 (P_{N-2} P_{N-3} ... P_2 P_1) Ĝ^(1)_11, indices decreasing left to right; N = 3: -Z^-1 P_1 Ĝ^(1)_11
        X = be.neg(_prod(be, [GNN] + [B.P(k) for k in range(N - 2, 0, -1)] + [B.pre(1, 1, 1)]))
    elif m == "M10":
        if lam_minus_1 is None:
            raise ValueError("M10 needs lam_minus_1 (Lambda - I)")
        Z = be.inv(GNN)
        X = be.inv(be.sub(Z, be.diag(lam_minus_1)))      # (Z - (Lambda - I))^-1
    elif m == "M12":
        X = G11                                           # algebraic identity
    elif m == "M13":
        X = G11
        for j in loop:
            X = be.add(X, be.scal(2, B.theta(j)))
    elif m == "M14":
        X = be.sub(G11, be.scal(2, B.pre(1, 1, 1)))
        sec = {"G_N1": be.neg(GN1), "G_1N": be.neg(G1N)}
    elif m == "M16":
        Zi3 = _prod(be, [GNN, GNN, GNN])                  # Z^-3 = Ĝ_NN^3
        X = be.add(B.pre(N - 1, 1, 1), _prod(be, [B.psi(N - 1), Zi3, B.chi(N - 1)]))
    else:  # pragma: no cover
        raise ValueError(m)
    return {"status": "computed", "designated": DESIGNATED[m], "X": X, "secondary": sec}


def reference_block(name, N, ny, G_full, backend):
    """The exact reference value of an endpoint block (G_11, G_N1, G_1N, G_NN) of the correct code."""
    a, b = {"G_11": (1, 1), "G_N1": (N, 1), "G_1N": (1, N), "G_NN": (N, N)}[name]
    return backend.block(G_full, a, b, ny)


def monomials(m, N):
    """Monomial structure of X_m - X_ref for Mag_m and B_tail (§4.3, §5.2): list of (|coefficient|, degree) of every
    monomial of the oracle formula of X_m *including* the reference block (coefficient 1, degree 1). Read off §4.2."""
    if not on_path(m, N) or m in ("M11", "M10"):
        raise OffPath(f"{m}: no polynomial monomial structure at N = {N}")
    L = N - 2                                             # number of loop slices
    ref = [(1, 1)]
    if m in ("M1", "M3", "M2"):
        return [(1, N)] + ref                             # Ĝ_NN Ĝ^(1)_11 P_2..P_{N-1}: N factors
    if m == "M4":
        return [(1, 1)] + [(1, 5), (1, 3)] * L + ref
    if m in ("M4L", "M4R"):
        return [(1, 1)] + [(1, 4), (1, 3)] * L + ref
    if m == "M5":
        return [(1, 1)] + [(1, 3), (1, 3)] * L + ref
    if m in ("M6", "M15", "M7"):
        return [(1, 2)] + ref
    if m == "M8":
        return [(1, 1), (1, 3), (1, 3)] + ref
    if m == "M9":
        return [(1, N)] + ref                             # Ĝ_NN P_{N-2}..P_1 Ĝ^(1)_11: N factors
    if m == "M12":
        return [(1, 1)] + ref
    if m == "M13":
        return [(1, 1)] + [(2, 3)] * L + ref
    if m == "M14":
        return [(1, 1), (2, 1)] + ref
    if m == "M16":
        return [(1, 1), (1, 5)] + ref                     # Ĝ^(N-1)_11 + psi Ĝ_NN^3 chi
    raise ValueError(m)


def monomial_factors(m, N, ny, G_full, G_prefix, backend):
    """The monomials of the oracle formula of X_m including the reference block (§4.3 Mag_m), as a list of
    (|coefficient|, [factor blocks]). Consistent with `monomials` (same count, coefficients and degrees).
    M10 (not a polynomial) is defined as [(1, [X_oracle]), (1, [Ĝ_NN])] — implementer choice C-7 (Mag for τ_abs,M10)."""
    be = backend
    if not on_path(m, N) or m == "M11":
        raise OffPath(f"{m} is off path at N = {N}")
    B = _Blocks(be, G_full, G_prefix, N, ny)
    GNN, G11 = B.full(N, N), B.full(1, 1)
    loop = list(range(2, N))
    ref = [(1, [reference_block(DESIGNATED[m], N, ny, G_full, be)])]
    if m in ("M1", "M3"):
        return [(1, [GNN, B.pre(1, 1, 1)] + [B.P(j) for j in loop])] + ref
    if m == "M2":
        return [(1, [B.P(j) for j in reversed(loop)] + [B.pre(1, 1, 1), GNN])] + ref
    if m == "M4":
        out = [(1, [G11])]
        for j in loop:
            P = B.P(j)
            out += [(1, [B.psi(j - 1), P, P, P, B.chi(j - 1)]), (1, [B.psi(j - 1), P, B.chi(j - 1)])]
        return out + ref
    if m in ("M4L", "M4R"):
        out = [(1, [G11])]
        for j in loop:
            P = B.P(j)
            out += [(1, [B.psi(j - 1), P, P, B.chi(j - 1)]), (1, [B.psi(j - 1), P, B.chi(j - 1)])]
        return out + ref
    if m == "M5":
        out = [(1, [G11])]
        for j in loop:
            out += [(1, [B.chi(j - 1), B.P(j), B.psi(j - 1)]), (1, [B.psi(j - 1), B.P(j), B.chi(j - 1)])]
        return out + ref
    if m in ("M6", "M15"):
        return [(1, [B.chi(N - 1), GNN])] + ref
    if m == "M7":
        return [(1, [GNN, B.psi(N - 1)])] + ref
    if m == "M8":
        return [(1, [G11]), (1, [B.chi(N - 1), GNN, B.psi(N - 1)]), (1, [B.psi(N - 1), GNN, B.chi(N - 1)])] + ref
    if m == "M9":
        return [(1, [GNN] + [B.P(k) for k in range(N - 2, 0, -1)] + [B.pre(1, 1, 1)])] + ref
    if m == "M10":
        return [(1, [None]), (1, [GNN])]                  # [None] = X_oracle, supplied by `magnitude`
    if m == "M12":
        return [(1, [G11])] + ref
    if m == "M13":
        return [(1, [G11])] + [(2, [B.psi(j - 1), B.P(j), B.chi(j - 1)]) for j in loop] + ref
    if m == "M14":
        return [(1, [G11]), (2, [B.pre(1, 1, 1)])] + ref
    if m == "M16":
        return [(1, [B.pre(N - 1, 1, 1)]), (1, [B.psi(N - 1), GNN, GNN, GNN, B.chi(N - 1)])] + ref
    raise ValueError(m)  # pragma: no cover


def magnitude(m, N, ny, G_full, G_prefix, backend, norm2, X_oracle=None):
    """Mag_m := Σ over the oracle monomials (incl. the reference block) of |coef|·Π‖factor‖₂ (§4.3).
    norm2: callable returning the operator 2-norm of a block. For M10, X_oracle (its oracle value) must be supplied."""
    total = 0
    for c, facs in monomial_factors(m, N, ny, G_full, G_prefix, backend):
        if m == "M10" and len(facs) == 1 and facs[0] is None:
            if X_oracle is None:
                raise ValueError("M10 magnitude needs X_oracle")
            facs = [X_oracle]
        p = c
        for f in facs:
            p = p * norm2(f)
        total = total + p
    return total
