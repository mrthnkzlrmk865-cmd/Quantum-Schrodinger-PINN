"""Stage 3 frozen-input preconditions (reference only): P5 N1-N5 and status list, mutation witnesses, targeted-lattice
lead-factor sensitivity, and the P4 coarse-enclosure omega >= 1. A failure here is PRECONDITION-UNMET."""
from fractions import Fraction as F

from cs2_t6_reference import rmat, rpoly, rlead, rcert
from . import inputs as I
from .common import Results, ref_blocks, t6_cases, recursion_diff, mids, yt6, eps
from .ychoices import Y_pert, Y_adv
from . import p5


def _swap_first(monos):
    (c, keys), rest = monos[0], monos[1:]
    k = list(keys)
    k[0], k[1] = k[1], k[0]
    return [(c, k)] + rest


def _swap_m5(monos):
    out = []
    for i, (c, keys) in enumerate(monos):
        k = list(keys)
        if i % 2 == 0:
            k[0], k[1] = k[1], k[0]
        out.append((c, k))
    return out


def run():
    R = Results("3 frozen-input preconditions")
    for it in p5.preconditions():
        R.add(*it)
    cases = list(t6_cases())
    # MU7 witnesses (per code path)
    for path, ms, sw in (("M1/M3", ("M1", "M3"), _swap_first), ("M2", ("M2",), _swap_first), ("M5", ("M5",), _swap_m5),
                         ("M6/M15", ("M6", "M15"), _swap_first), ("M7", ("M7",), _swap_first), ("M8", ("M8",), _swap_first)):
        ok = any(rpoly.evaluate(sw(rpoly.monomials(m, N)), ref_blocks(n, N, s)) != rpoly.F_ref(m, N, ref_blocks(n, N, s))
                 for n, N, s, m in cases if m in ms)
        R.add(f"MU7 witness {path}: swapped monomial changes F on some T6.1 input", ok, "PRECONDITION-UNMET")
    for path, ms, key in (("M1/M3", ("M1", "M3"), "G_N1"), ("M2", ("M2",), "G_1N"), ("M6/M15", ("M6", "M15"), "G_N1"),
                          ("M7", ("M7",), "G_1N")):
        ok = any(not rmat.is_zero(ref_blocks(n, N, s)[key]) for n, N, s, m in cases if m in ms)
        R.add(f"MU8 witness {path}: single block {key} != 0 on some input", ok, "PRECONDITION-UNMET")
    ok = any(not rmat.is_zero(rpoly.evaluate(rpoly.monomials("M5", 4)[2:], ref_blocks(n, 4, s))) for n, N, s, m in cases
             if m == "M5" and N == 4)
    R.add("MU9 witness: the j = 3 term of M5 at N = 4 is nonzero on some input", ok, "PRECONDITION-UNMET")
    ok = False
    for name, N in (("X1", 3), ("X2", 4)):
        Ys = Y_pert(name, N, F(1, 2))
        Bt, _ = rcert.blocks_from_Y(Ys, N, len(I.LATTICES[name]["gammas"]))
        B = ref_blocks(name, N, F(1, 2))
        ok |= any(rpoly.F_ref(m, N, B) != rpoly.F_ref(m, N, Bt) for m in I.ORDER_MUTANTS if I.on_path(m, N))
    R.add("MU10 witness: F(B) - F(B~) != 0 for Y_pert on some T6.2 input", ok, "PRECONDITION-UNMET")
    Ya = Y_adv()
    Bt, _ = rcert.blocks_from_Y(Ya, 2, 2)
    R.add("MU11 witness: F_M6(B~_adv) = 0 exactly and D*_true > 0 (X1 N=2 s=2^-19)",
          rmat.is_zero(rpoly.F_ref("M6", 2, Bt)) and rmat.fro2(recursion_diff("X1", 2, F(1, 2**19), "M6")) > 0,
          "PRECONDITION-UNMET")
    # targeted gamma = 5/4 lattices: ||A_e||_F^2 = k u^2 exactly, correct eps passes, MU1-MU5 formulas fail (§2.4.3.9)
    u = I.T * F(1, 2**257)
    for name, N, kinds in (("G1", 1, {"full": 4}), ("G2", 2, {"full": 4, ("prefix", 1): 2}),
                           ("G3", 3, {"full": 2, ("prefix", 1): 1, ("prefix", 2): 1})):
        At = rcert.true_matrices(I.LATTICES[name], N, F(1))
        Am = mids(name, N, F(1))
        for key, k in kinds.items():
            ae2 = rmat.fro2(rmat.sub(At[key], Am[key]))
            e = eps(name, N, "full" if key == "full" else "prefix")
            R.add(f"targeted {name} {key}: ||A_e||_F^2 = {k} u^2 exactly and correct eps passes",
                  ae2 == k * u * u and e * e >= ae2, "PRECONDITION-UNMET")
    R.add("targeted MU1-MU5 formulas fail on G1-G3 (exact): ub(2u^2)^2 < 4u^2, ub(u^2)^2 < 2u^2, lb(2u^2)^2 < 2u^2",
          rlead.ub(2 * u * u) ** 2 < 4 * u * u and rlead.ub(u * u) ** 2 < 2 * u * u and rlead.lb(2 * u * u) ** 2 < 2 * u * u,
          "PRECONDITION-UNMET")
    # P4: the reference shows omega(full) >= 1 for C8 at ell = 4 with Y_T6 (and every point is the same matrix)
    s = F(1, 2)
    A4 = mids("C8", 2, s, 4)["full"]
    c = rcert.cert_point_ref(A4, yt6("C8", 2, s, 4)["full"], eps("C8", 2, "full", 4))
    R.add("P4 reference omega(full) >= 1 at ell = 4 (and > 1.39 per §2.4.3.9)", c["omega"] >= 1 and c["omega"] > F(139, 100),
          "PRECONDITION-UNMET", f"omega = {float(c['omega']):.4f}")
    return R
