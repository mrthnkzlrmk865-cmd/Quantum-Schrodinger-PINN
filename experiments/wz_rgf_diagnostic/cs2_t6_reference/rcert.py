"""Independent implementation of the §2.4.1/§5.6 certificate quantities (eps_A, omega, a, h, eta, beta, E_m, U_m, G_lo)
from reference modules only (own ub/lb from rlead). Used for P5 (N4, N5, value agreement) and for mutation witnesses."""
from fractions import Fraction
from .qexact import RC
from . import rmat, rinv, rlead, rlattice, rpoly

T = rlattice.T


def enclosures(lat, ell=rlead.ELL):
    return [rlead.enclosure(g, ell) for g in lat["gammas"]]


def eps_A_ref(lat, N, kind, ell=rlead.ELL):
    te = [T * enc["e"] for enc in enclosures(lat, ell)]
    if kind == "prefix":
        return rlead.ub(sum((x * x for x in te), Fraction(0)))
    if N == 1:
        return rlead.ub(sum(((2 * x) * (2 * x) for x in te), Fraction(0)))
    return rlead.ub(2 * sum((x * x for x in te), Fraction(0)))


def matrix_keys(N):
    return ["full"] + [("prefix", j) for j in range(1, N)]


def mid_matrices(lat, N, s, ell=rlead.ELL):
    lam_mid = [enc["lam_mid"] for enc in enclosures(lat, ell)]
    out = {"full": rlattice.M_phys(lat, N, s, "full", lam=lam_mid)}
    for j in range(1, N):
        out[("prefix", j)] = rlattice.M_phys(lat, N, s, "prefix", j=j, lam=lam_mid)
    return out


def true_matrices(lat, N, s):
    lam = [rlead.exact_root(g) for g in lat["gammas"]]
    out = {"full": rlattice.M_phys(lat, N, s, "full", lam=lam)}
    for j in range(1, N):
        out[("prefix", j)] = rlattice.M_phys(lat, N, s, "prefix", j=j, lam=lam)
    return out


def y_t6(A_mid):
    return rlead.round_matrix(rinv.inverse(A_mid), rlead.ELL_Y)


def cert_point_ref(A_mid, Y, eps):
    n = len(A_mid)
    Xi = rmat.sub(rmat.eye(n), rmat.matmul(A_mid, Y))
    f_Xi, f_Y, f_mid = rmat.fro2(Xi), rmat.fro2(Y), rmat.fro2(A_mid)
    omega = rlead.ub(f_Xi) + eps * rlead.ub(f_Y)
    h = rlead.ub(f_mid) + eps
    out = dict(f_Xi=f_Xi, f_Y=f_Y, f_mid=f_mid, omega=omega, h=h, certified=omega < 1)
    if omega < 1:
        out["a"] = rlead.ub(f_Y) / (1 - omega)
        out["eta"] = T * rlead.ub(f_Y) * omega / (1 - omega)
    return out


def required_prefixes(m, N):
    if m in ("M1", "M2", "M3", "M5"):
        return list(range(1, N))
    return [N - 1]


def blocks_from_Y(Ys, N, ny):
    """Blocks of t*Y: full -> G_NN, G_N1, G_1N; prefix j -> P_j = -(tY)_jj, chi_j, psi_j; G1_11 from prefix 1."""
    out, owner = {}, {}
    tY = rmat.scale(RC(T), Ys["full"])
    for k, (a, b) in (("G_NN", (N, N)), ("G_N1", (N, 1)), ("G_1N", (1, N))):
        out[k] = rmat.block(tY, a, b, ny)
        owner[k] = "full"
    for j in range(1, N):
        if ("prefix", j) not in Ys:
            continue
        tYp = rmat.scale(RC(T), Ys[("prefix", j)])
        out[("P", j)] = rmat.neg(rmat.block(tYp, j, j, ny))
        out[("chi", j)] = rmat.block(tYp, j, 1, ny)
        out[("psi", j)] = rmat.block(tYp, 1, j, ny)
        owner[("P", j)] = owner[("chi", j)] = owner[("psi", j)] = ("prefix", j)
        if j == 1:
            out["G1_11"] = rmat.block(tYp, 1, 1, ny)
            owner["G1_11"] = ("prefix", 1)
    return out, owner


def pointwise_ref(m, N, ny, refs, monos=None):
    """refs: key -> (A_mid, Y, eps). Returns the §5.6 quantities; raises if a required matrix is uncertified."""
    certs = {k: cert_point_ref(*refs[k]) for k in ["full"] + [("prefix", j) for j in required_prefixes(m, N)]}
    for k, c in certs.items():
        if not c["certified"]:
            raise ArithmeticError(f"reference matrix {k} uncertified (omega >= 1)")
    Ys = {k: refs[k][1] for k in certs}
    B, owner = blocks_from_Y(Ys, N, ny)
    monos = rpoly.monomials(m, N) if monos is None else monos
    beta = {}
    E_m = Fraction(0)
    for c, keys in monos:
        prod_up, prod = Fraction(1), Fraction(1)
        for k in keys:
            if k not in beta:
                beta[k] = rlead.ub(rmat.fro2(B[k]))
            prod_up *= beta[k] + certs[owner[k]]["eta"]
            prod *= beta[k]
        E_m += abs(c) * (prod_up - prod)
    F = rpoly.evaluate(monos, B)
    U = rlead.ub(rmat.fro2(F)) + rlead.ub(Fraction(ny)) * E_m
    full = certs["full"]
    h_W = full["h"]
    G_lo = max(T / h_W, T * rlead.lb(rmat.max_abs2(refs["full"][1])) - full["eta"])
    return dict(U=U, G_lo=G_lo, E_m=E_m, beta=beta, eta={k: c["eta"] for k, c in certs.items()}, certs=certs, F=F,
                blocks=B)
