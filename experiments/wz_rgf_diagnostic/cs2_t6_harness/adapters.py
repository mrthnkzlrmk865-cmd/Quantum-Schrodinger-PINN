"""Boundary B1 (reference RC -> production QC pairs) and B2 (production output -> harness, exact types only)."""
from fractions import Fraction

from cs2_t6_reference import rmat, rlattice
from cs2_t6_reference.qexact import RC
from . import inputs as I
from .common import mids, eps


class InterfaceFail(TypeError):
    pass


def exact_out(x, what="output"):
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise InterfaceFail(f"{what}: production returned {type(x).__name__}, exact int/Fraction required")
    return Fraction(x)


def to_prod(X):
    return tuple(tuple((z.re, z.im) for z in r) for r in X)


def from_prod(Q, what="matrix"):
    try:
        return tuple(tuple(RC(exact_out(z[0], what), exact_out(z[1], what)) for z in r) for r in Q)
    except (TypeError, IndexError) as exc:
        raise InterfaceFail(f"{what}: {exc}")


def delta_prod(name, N, s):
    return [to_prod(D) for D in rlattice.delta_blocks(I.LATTICES[name], N, s)]


def prod_mid(prod, name, N, s, key, ell=256):
    lat = I.LATTICES[name]
    kind, j = ("full", None) if key == "full" else ("prefix", key[1])
    return prod.assemble_A_mid(delta_prod(name, N, s), lat["gammas"], kind, j, ell)


def prod_eps(prod, name, N, key, ell=256):
    lat = I.LATTICES[name]
    return exact_out(prod.eps_A(lat["gammas"], N, "full" if key == "full" else "prefix", ell), "eps_A")


def check_mid(R, prod, name, N, s, ell=256):
    """Production A_mid (Delta form) must equal the reference A_mid (H form) exactly; returns production refs base."""
    out = {}
    ok_all = True
    for key, Aref in mids(name, N, s, ell).items():
        Ap = from_prod(prod_mid(prod, name, N, s, key, ell), "A_mid")
        ok = Ap == Aref
        ok_all &= ok
        e_p, e_r = prod_eps(prod, name, N, key, ell), eps(name, N, "full" if key == "full" else "prefix", ell)
        out[key] = (to_prod(Ap), e_p)
    R.add(f"A_mid production == reference [{name} N={N} s={s} ell={ell}]", ok_all, "CERT-VIOLATION")
    return out


def prod_refs(base, Ys, keys=None):
    keys = list(base) if keys is None else keys
    return {k: (base[k][0], to_prod(Ys[k]), base[k][1]) for k in keys if k in Ys}
