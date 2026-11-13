"""T6.2: eta_A >= t*2^-40 (Y_pert) and ||F(B) - F(B~)||_2 <= E_m (Y_pert and Y_exact)."""
from fractions import Fraction as F
from cs2_t6_reference import rmat, rpsd, rpoly, rcert
from . import inputs as I
from .common import Results, ref_blocks, ref_certified, required_keys, expect_status
from .adapters import check_mid, prod_refs, exact_out
from .ychoices import Y_pert, Y_exact

CASES = (("X1", 3), ("X2", 4))


def run(prod, R=None):
    R = R or Results("T6.2")
    s = F(1, 2)
    for name, N in CASES:
        base = check_mid(R, prod, name, N, s)
        ny = len(I.LATTICES[name]["gammas"])
        B = ref_blocks(name, N, s)
        for yname, Ys in (("Y_pert", Y_pert(name, N, s)), ("Y_exact", Y_exact(name, N, s))):
            Bt, _ = rcert.blocks_from_Y(Ys, N, ny)
            for m in I.ORDER_MUTANTS:
                if not I.on_path(m, N):
                    continue
                res = prod.pointwise_certificate(m, N, ny, prod_refs(base, Ys))
                tag = f"[{name} N={N} {yname} {m}]"
                if not expect_status(R, res, ref_certified(name, N, s, required_keys(m, N), Ys), f"T6.2 {tag}"):
                    continue
                if yname == "Y_pert":
                    R.add(f"T6.2 eta >= t*2^-40 {tag}",
                          all(exact_out(v, "eta") >= I.T * F(1, 2**40) for v in res["eta"].values()), "CERT-VIOLATION")
                D = rmat.sub(rpoly.F_ref(m, N, B), rpoly.F_ref(m, N, Bt))
                R.add(f"T6.2 ||F(B)-F(B~)||_2 <= E_m {tag}", rpsd.norm2_le(D, exact_out(res["E_m"], "E_m")),
                      "CERT-VIOLATION")
    return R
