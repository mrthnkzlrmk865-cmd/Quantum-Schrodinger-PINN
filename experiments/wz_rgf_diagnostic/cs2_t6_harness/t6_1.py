"""T6.1: F_poly on exact blocks == exact mutant recursion designated block - exact reference block."""
from cs2_t6_reference import rmat
from . import inputs as I
from .common import Results, ref_blocks, recursion_diff, t6_cases
from .adapters import to_prod, from_prod

KEYS = ("G_NN", "G_N1", "G_1N", "G1_11")


def prod_blocks(b, N):
    out = {k: to_prod(b[k]) for k in KEYS if k in b}
    for j in range(1, N):
        for t in ("P", "chi", "psi"):
            out[(t, j)] = to_prod(b[(t, j)])
    return out


def run(prod, R=None):
    R = R or Results("T6.1")
    for name, N, s, m in t6_cases():
        b = ref_blocks(name, N, s)
        got = from_prod(prod.F_poly(m, N, prod_blocks(b, N)), "F_poly")
        R.add(f"T6.1 [{name} N={N} s={s} {m}]", got == recursion_diff(name, N, s, m), "CERT-VIOLATION")
    return R
