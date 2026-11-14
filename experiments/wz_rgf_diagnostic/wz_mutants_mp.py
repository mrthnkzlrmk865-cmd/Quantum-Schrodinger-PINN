"""Arbitrary-precision (mpmath) W/Z endpoint blocks with code mutants.

Clean-room implementation written ONLY from ``wz_mutants_mp_MATERIALS.md``.
Implements ``mp_endpoint_blocks(Delta, lam_minus_1, mutant=None, dps=60)``.
"""

import mpmath

PROVENANCE = (
    "wz_mutants_mp.py was written in a clean-room context solely from the materials file "
    "experiments/wz_rgf_diagnostic/wz_mutants_mp_MATERIALS.md (correct-code listing, change column, "
    "N = 1 rules and interface); no other repository source was read."
)

# Section-3 labels, in the order in which they appear in Section 3.
MUTANTS = ("M1", "M2", "M3", "M4", "M4L", "M4R", "M5", "M6", "M7", "M8", "M9", "M10",
           "M12", "M13", "M14", "M15", "M16", "M11")


def _to_matrix(ctx, x):
    if hasattr(x, "tolist"):
        x = x.tolist()
    rows = [[ctx.mpc(v) for v in row] for row in x]
    return ctx.matrix(rows)


def _solve(ctx, X, B):
    """X^-1 B, column by column with ctx.lu_solve."""
    n = X.rows
    out = ctx.matrix(n, B.cols)
    for c in range(B.cols):
        col = ctx.matrix([B[r, c] for r in range(B.rows)])
        x = ctx.lu_solve(X, col)
        for r in range(n):
            out[r, c] = x[r]
    return out


def mp_endpoint_blocks(Delta, lam_minus_1, mutant=None, dps=60):
    if mutant is not None and mutant not in MUTANTS:
        raise ValueError("unknown mutant label: %r" % (mutant,))
    m = mutant
    ctx = mpmath.MPContext()
    ctx.dps = dps

    D = [_to_matrix(ctx, d) for d in Delta]
    N = len(D)
    if N < 1:
        raise ValueError("Delta must contain at least one slice")
    lm = [ctx.mpc(v) for v in lam_minus_1]
    ny = len(lm)
    I = ctx.eye(ny)
    Lm = ctx.matrix(ny, ny)
    for k in range(ny):
        Lm[k, k] = lm[k]
    for d in D:
        if d.rows != ny or d.cols != ny:
            raise ValueError("Delta blocks must be N_y x N_y with N_y = len(lam_minus_1)")

    trace = []

    def rec(name, M):
        trace.append((name, M.copy()))

    def inv_I_minus_W_times_W(W):
        # (I - W)^-1 W
        if m == "M12":
            return _solve(ctx, (I - W).T, W.T).T
        return _solve(ctx, I - W, W)

    if N == 1:
        if m == "M11":
            Z = 2 * D[0] + Lm
        else:
            Z = 2 * D[0] + 2 * Lm
        rec("Z", Z)
        Zi = _solve(ctx, Z, I)
        rec("Zi", Zi)
        G11 = Zi.copy(); GN1 = Zi.copy(); G1N = Zi.copy(); GNN = Zi.copy()
        rec("G11", G11); rec("GN1", GN1); rec("G1N", G1N); rec("GNN", GNN)
        return {"Ghat11": G11, "GhatN1": GN1, "Ghat1N": G1N, "GhatNN": GNN,
                "trace": trace, "Z": Z, "Zi": Zi}

    # init
    W = 2 * D[0] + Lm
    rec("W_1", W)
    if m == "M14":
        G11 = _solve(ctx, I - W, I)
    else:
        G11 = _solve(ctx, W - I, I)
    rec("G11^(1)", G11)
    Gj1 = G11.copy()
    G1j = G11.copy()

    # loop over slices s = 2..N-1
    for s in range(2, N):
        Ds = D[s - 1]
        if m == "M9":
            P = _solve(ctx, I - W, I)  # from W_old, before the W update
            rec("P_%d" % s, P)
            W = 2 * Ds + inv_I_minus_W_times_W(W)
            rec("W_%d" % s, W)
        else:
            W = 2 * Ds + inv_I_minus_W_times_W(W)
            rec("W_%d" % s, W)
            P = _solve(ctx, I - W, I)
            rec("P_%d" % s, P)

        def upd_G11(G11, Gj1, G1j):
            if m == "M5":
                return G11 - Gj1 * P * G1j
            if m == "M13":
                return G11 + G1j * P * Gj1
            return G11 - G1j * P * Gj1

        def upd_Gj1(Gj1):
            if m in ("M1", "M3"):
                return Gj1 * P
            return P * Gj1

        def upd_G1j(G1j):
            if m in ("M2", "M3"):
                return P * G1j
            return G1j * P

        if m == "M4":
            Gj1 = upd_Gj1(Gj1); rec("G_%d1^(%d)" % (s, s), Gj1)
            G1j = upd_G1j(G1j); rec("G_1%d^(%d)" % (s, s), G1j)
            G11 = upd_G11(G11, Gj1, G1j); rec("G11^(%d)" % s, G11)
        elif m == "M4L":
            Gj1 = upd_Gj1(Gj1); rec("G_%d1^(%d)" % (s, s), Gj1)
            G11 = upd_G11(G11, Gj1, G1j); rec("G11^(%d)" % s, G11)
            G1j = upd_G1j(G1j); rec("G_1%d^(%d)" % (s, s), G1j)
        elif m == "M4R":
            G1j = upd_G1j(G1j); rec("G_1%d^(%d)" % (s, s), G1j)
            G11 = upd_G11(G11, Gj1, G1j); rec("G11^(%d)" % s, G11)
            Gj1 = upd_Gj1(Gj1); rec("G_%d1^(%d)" % (s, s), Gj1)
        else:
            G11 = upd_G11(G11, Gj1, G1j); rec("G11^(%d)" % s, G11)
            Gj1 = upd_Gj1(Gj1); rec("G_%d1^(%d)" % (s, s), Gj1)
            G1j = upd_G1j(G1j); rec("G_1%d^(%d)" % (s, s), G1j)

    # closing
    if m == "M10":
        Z = 2 * D[N - 1] + inv_I_minus_W_times_W(W)
    else:
        Z = 2 * D[N - 1] + Lm + inv_I_minus_W_times_W(W)
    rec("Z", Z)
    Zi = _solve(ctx, Z, I)
    rec("Zi", Zi)

    def mk_GN1():
        if m in ("M6", "M15"):
            return -(Gj1 * Zi)
        return -(Zi * Gj1)

    def mk_G1N():
        if m in ("M7", "M15"):
            return -(Zi * G1j)
        return -(G1j * Zi)

    if m == "M16":
        GN1 = mk_GN1(); rec("GN1", GN1)
        G1N = mk_G1N(); rec("G1N", G1N)
        G11 = G11 + G1N * Zi * GN1; rec("G11", G11)
    else:
        if m == "M8":
            G11 = G11 + Gj1 * Zi * G1j
        else:
            G11 = G11 + G1j * Zi * Gj1
        rec("G11", G11)
        GN1 = mk_GN1(); rec("GN1", GN1)
        G1N = mk_G1N(); rec("G1N", G1N)
    GNN = Zi.copy()
    rec("GNN", GNN)

    return {"Ghat11": G11, "GhatN1": GN1, "Ghat1N": G1N, "GhatNN": GNN,
            "trace": trace, "Z": Z, "Zi": Zi}
