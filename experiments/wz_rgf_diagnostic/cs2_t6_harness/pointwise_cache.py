"""Production pointwise results for the T6.5/T6.6 cases (Y_T6), computed once per production module object."""
from . import inputs as I
from .common import t6_cases
from .adapters import check_mid, prod_refs
from .ychoices import Y_T6


def compute(R, prod):
    out = {}
    bases = {}
    for name, N, s, m in t6_cases():
        if (name, N, s) not in bases:
            bases[(name, N, s)] = check_mid(R, prod, name, N, s)
        ny = len(I.LATTICES[name]["gammas"])
        out[(name, N, s, m)] = prod.pointwise_certificate(m, N, ny, prod_refs(bases[(name, N, s)], Y_T6(name, N, s)))
    return out
