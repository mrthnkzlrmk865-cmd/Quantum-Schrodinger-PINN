"""Independent references for T1-T5 (written here; no production helper is imported or called):
a Fraction-only isqrt enclosure, a Fraction-only second assembly routine from M = E - H - Sigma_L (+) Sigma_R, and exact
closed-form singular-value checks."""
from fractions import Fraction as F
from math import isqrt

T = F(32)
DX2 = F(1, 64)
ELL = 256


def enclosure(gamma, ell=ELL):
    g = F(gamma)
    if abs(g) < 1:
        br, xi = "open", 1 - g * g
    elif g > 1:
        br, xi = "closed", g * g - 1
    else:
        raise ValueError("excluded gamma")
    a = xi * 4 ** ell
    fl = a.numerator // a.denominator
    ce = -((-a.numerator) // a.denominator)
    lo = F(isqrt(fl), 2 ** ell)
    hi = F(isqrt(ce) + 1, 2 ** ell)
    mu = (lo + hi) / 2
    return dict(branch=br, xi=xi, lo=lo, hi=hi, mu=mu, e=(hi - lo) / 2,
                lam=(g - mu, F(0)) if br == "closed" else (g, mu))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def neg(a):
    return (-a[0], -a[1])


def _blank(n):
    return [[(F(0), F(0)) for _ in range(n)] for _ in range(n)]


def _finish(M):
    return tuple(tuple(r) for r in M)


def second_assembly_L(E, E_perp, V1, s, n, full):
    """M_ab = E delta_ab - H_ab - Sigma_ab with H = blocktridiag(2t I + diag(E_perp) + s V_j, -t I), Sigma = -t lambda~."""
    ny = len(E_perp)
    M = _blank(n * ny)
    for j in range(n):
        for a in range(ny):
            for b in range(ny):
                H = (s * V1[j][a][b][0], s * V1[j][a][b][1])
                if a == b:
                    H = (H[0] + 2 * T + E_perp[a], H[1])
                    M[j * ny + a][j * ny + b] = (E - H[0], -H[1])
                else:
                    M[j * ny + a][j * ny + b] = neg(H)
        if j + 1 < n:
            for a in range(ny):
                M[j * ny + a][(j + 1) * ny + a] = (T, F(0))       # -H_{j,j+1} = +t
                M[(j + 1) * ny + a][j * ny + a] = (T, F(0))
    sig = [neg(cmul((T, F(0)), enclosure(1 - (E - e) * DX2)["lam"])) for e in E_perp]   # Sigma = -t lambda~
    for a in range(ny):
        M[a][a] = cadd(M[a][a], neg(sig[a]))
        if full:
            k = (n - 1) * ny + a
            M[k][k] = cadd(M[k][k], neg(sig[a]))
    return _finish(M)


def second_assembly_R_W(Delta, dlead, n, full):
    """M/t has diagonal blocks 2 Delta_j - 2I (+ Lambda~ on the lead slices) and off-diagonal I; M = t (M/t)."""
    ny = len(dlead)
    M = _blank(n * ny)
    for j in range(n):
        for a in range(ny):
            for b in range(ny):
                d = Delta[j][a][b]
                v = (2 * d[0] - (2 if a == b else 0), 2 * d[1])
                M[j * ny + a][j * ny + b] = (T * v[0], T * v[1])
        if j + 1 < n:
            for a in range(ny):
                M[j * ny + a][(j + 1) * ny + a] = (T, F(0))
                M[(j + 1) * ny + a][j * ny + a] = (T, F(0))
    lam = [enclosure(1 - d)["lam"] for d in dlead]
    for a in range(ny):
        M[a][a] = cadd(M[a][a], cmul((T, F(0)), lam[a]))
        if full:
            k = (n - 1) * ny + a
            M[k][k] = cadd(M[k][k], cmul((T, F(0)), lam[a]))
    return _finish(M)


def second_assembly_R_orig(A_orig, sigma, n, full):
    ny = len(sigma)
    M = _blank(n * ny)
    for j in range(n):
        for a in range(ny):
            for b in range(ny):
                M[j * ny + a][j * ny + b] = A_orig[j][a][b]
        if j + 1 < n:
            for a in range(ny):
                M[j * ny + a][(j + 1) * ny + a] = (T, F(0))
                M[(j + 1) * ny + a][j * ny + a] = (T, F(0))
    for a in range(ny):
        M[a][a] = cadd(M[a][a], neg(sigma[a]))
        if full:
            k = (n - 1) * ny + a
            M[k][k] = cadd(M[k][k], neg(sigma[a]))
    return _finish(M)


def sigma_max_2x2_le(X, a):
    """Exact a >= sigma_max(X) for a 2x2 complex X (§2.4.2 T1): with F = ||X||_F^2 and D = |det X|^2,
    sigma_max^2 = (F + sqrt(F^2 - 4D))/2, so a >= sigma_max iff 2a^2 - F >= 0 and (2a^2 - F)^2 >= F^2 - 4D."""
    (p, q), (r, s) = X
    Fn = sum(z[0] ** 2 + z[1] ** 2 for z in (p, q, r, s))
    det = (p[0] * s[0] - p[1] * s[1] - (q[0] * r[0] - q[1] * r[1]), p[0] * s[1] + p[1] * s[0] - (q[0] * r[1] + q[1] * r[0]))
    D = det[0] ** 2 + det[1] ** 2
    a = F(a)
    return 2 * a * a - Fn >= 0 and (2 * a * a - Fn) ** 2 >= Fn * Fn - 4 * D
