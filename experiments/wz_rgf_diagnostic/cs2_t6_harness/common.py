"""Shared harness utilities: result records and cached reference computations."""
from fractions import Fraction
from functools import lru_cache

from cs2_t6_reference import rmat, rlead, rlattice, rblocks, rmutants, rpoly, rcert
from cs2_t6_reference.qexact import RC
from . import inputs as I


class Results:
    def __init__(self, stage):
        self.stage = stage
        self.items = []

    def add(self, name, ok, cls=None, detail=""):
        self.items.append(dict(name=name, ok=bool(ok), cls=None if ok else cls, detail=str(detail)))
        return ok

    @property
    def failures(self):
        return [x for x in self.items if not x["ok"]]

    def summary(self):
        return dict(stage=self.stage, total=len(self.items), failed=len(self.failures))


def true_lams(name):
    return [rlead.exact_root(g) for g in I.LATTICES[name]["gammas"]]


@lru_cache(maxsize=None)
def ref_blocks(name, N, s):
    return rblocks.blocks(I.LATTICES[name], N, s, true_lams(name))


@lru_cache(maxsize=None)
def recursion(name, N, s, m):
    lat = I.LATTICES[name]
    D = rlattice.delta_blocks(lat, N, s)
    return rmutants.endpoint_blocks(D, [l - RC(1) for l in true_lams(name)], m)


@lru_cache(maxsize=None)
def mids(name, N, s, ell=rlead.ELL):
    return rcert.mid_matrices(I.LATTICES[name], N, s, ell)


@lru_cache(maxsize=None)
def ytrue(name, N, s):
    b = ref_blocks(name, N, s)
    out = {"full": b["Ainv"]}
    for j in range(1, N):
        out[("prefix", j)] = b["prefix"][j]["Ainv"]
    return out


@lru_cache(maxsize=None)
def yt6(name, N, s, ell=rlead.ELL):
    return {k: rcert.y_t6(A) for k, A in mids(name, N, s, ell).items()}


def recursion_diff(name, N, s, m):
    rec = recursion(name, N, s, m)
    key = {"G_N1": "GN1", "G_1N": "G1N", "G_11": "G11"}[rpoly.DESIGNATED[m]]
    return rmat.sub(rec[key], ref_blocks(name, N, s)[rpoly.DESIGNATED[m]])


def eps(name, N, kind, ell=rlead.ELL):
    return rcert.eps_A_ref(I.LATTICES[name], N, kind, ell)


def t6_cases():
    """(lattice, N, s, mutant) of T6.1/T6.5."""
    for name in I.MAIN:
        for N in I.NS[name]:
            for s in I.S6:
                for m in I.ORDER_MUTANTS:
                    if I.on_path(m, N):
                        yield name, N, s, m


def all_cert_cases():
    """(lattice, N, s) of T6.3/T6.4."""
    for name in I.MAIN + I.TARGETED:
        for N in I.NS[name]:
            for s in I.S6:
                yield name, N, s


def ref_certified(name, N, s, keys, Ys, ell=256):
    """Reference-side expectation: are all listed matrices certified (omega < 1) with these Y? (reference A_mid, eps)."""
    A = mids(name, N, s, ell)
    for k in keys:
        kind = "full" if k == "full" else "prefix"
        if not rcert.cert_point_ref(A[k], Ys[k], eps(name, N, kind, ell))["certified"]:
            return False
    return True


def required_keys(m, N):
    return ["full"] + [("prefix", j) for j in rcert.required_prefixes(m, N)]


def expect_status(R, res, ref_ok, tag, accept_bounds=False):
    """Production status 'computed' is expected iff the reference certifies every required matrix.
    Returns True if the bound fields may be checked."""
    if res.get("status") == "computed":
        return True
    if accept_bounds and res.get("reason") in ("invalid_lower_bound", "invalid_upper_bound"):
        return True  # bounds present; the bound checks classify them (frozen-input invariant)
    if ref_ok:
        R.add(f"computed {tag}", False, "CERT-VIOLATION", f"production {res.get('status')}/{res.get('reason')} "
              "although the reference certifies every required matrix")
    else:
        R.add(f"computed {tag}", False, "PRECONDITION-UNMET", "reference does not certify a required matrix")
    return False
