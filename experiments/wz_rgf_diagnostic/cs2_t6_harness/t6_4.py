"""T6.4 exact: A_e zero off the leads, eps_A^2 >= ||A_e||_F^2, ||I - A_true Y||_2 <= omega (Y_T6); includes G1-G3."""
from cs2_t6_reference import rmat, rpsd, rcert
from . import inputs as I
from .common import Results, all_cert_cases
from .adapters import check_mid, exact_out, to_prod
from .ychoices import Y_T6


def run(prod, R=None, cases=None):
    R = R or Results("T6.4")
    for name, N, s in (cases or list(all_cert_cases())):
        base = check_mid(R, prod, name, N, s)
        ny = len(I.LATTICES[name]["gammas"])
        At = rcert.true_matrices(I.LATTICES[name], N, s)
        Ys = Y_T6(name, N, s)
        for key, (Ap, e) in base.items():
            from .adapters import from_prod
            Ae = rmat.sub(At[key], from_prod(Ap))
            n = len(Ae)
            last = N if key == "full" else key[1]
            lead = {(p, p) for p in range(ny)} | ({((last - 1) * ny + p,) * 2 for p in range(ny)} if key == "full" else set())
            t = f"[{name} N={N} s={s} {key}]"
            R.add(f"T6.4 A_e zero off leads {t}",
                  all(Ae[i][j].is_zero() for i in range(n) for j in range(n) if (i, j) not in lead), "CERT-VIOLATION")
            e = exact_out(e, "eps_A")
            R.add(f"T6.4 eps_A^2 >= ||A_e||_F^2 {t}", e * e >= rmat.fro2(Ae), "CERT-VIOLATION")
            c = prod.certify_point(Ap, to_prod(Ys[key]), e)
            Rm = rmat.sub(rmat.eye(n), rmat.matmul(At[key], Ys[key]))
            R.add(f"T6.4 ||I - A_true Y||_2 <= omega {t}", rpsd.norm2_le(Rm, exact_out(c["omega"], "omega")),
                  "CERT-VIOLATION")
    return R
