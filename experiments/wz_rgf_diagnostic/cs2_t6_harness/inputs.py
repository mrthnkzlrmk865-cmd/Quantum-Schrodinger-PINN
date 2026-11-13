"""Literal frozen T6 inputs of §2.4.3.1, §2.4.3.2, §2.4.3.6 and §2.4.3.8 (Revision 3.7). Exact rationals only."""
from fractions import Fraction as F

T = F(32)
DX2 = F(1, 64)
ELL = 256
ELL_Y = 512
THRESHOLD_LITERAL = F(1, 10**10)

g12, g13, g25, gP, gE = F(5, 4), F(5, 3), F(29, 20), F(3, 5), F(1025, 1024)

VS = [[[F(3)]], [[F(-5, 2)]], [[F(1)]], [[F(-2)]]]
VR = [[[F(2), F(3)], [F(3), F(-1)]],
      [[F(-4), F(1)], [F(1), F(5, 2)]],
      [[F(1, 2), F(-2)], [F(-2), F(3)]],
      [[F(-3), F(1, 2)], [F(1, 2), F(-1, 2)]]]
VH = [[[(F(1), F(0)), (F(2), F(1))], [(F(2), F(-1)), (F(-2), F(0))]],
      [[(F(3), F(0)), (F(-1), F(2))], [(F(-1), F(-2)), (F(1, 2), F(0))]],
      [[(F(-2), F(0)), (F(1, 2), F(-3, 2))], [(F(1, 2), F(3, 2)), (F(1), F(0))]],
      [[(F(1, 2), F(0)), (F(-3), F(1, 2))], [(F(-3), F(-1, 2)), (F(-1), F(0))]]]
VD = [[[v[0][0], F(0)], [F(0), v[1][1]]] for v in VR]

S6 = [F(1), F(1, 2), F(1, 2**19), F(1, 2**40)]
S_TEST = sorted({F(1, 2**k) for k in range(19, 41)} | {F(3, 2**k) for k in range(21, 42)}, reverse=True)
assert len(S_TEST) == 43


def lat(gammas, V):
    return {"gammas": list(gammas), "V": V}


LATTICES = {
    "X0": lat([gP], VS), "X1": lat([g13, gP], VR), "X2": lat([g25, gP], VH), "XD": lat([g12, gP], VD),
    "G1": lat([g12], VS[:1]), "G2": lat([g12, g12], VR[:2]), "G3": lat([g12], VS[:3]),
    "X1n": lat([g13, gP], [[[x / 128 for x in row] for row in VR[j]] for j in range(2)]),
    "C8": lat([gE], [[[F(0)]], [[F(0)]]]),
}
NS = {"X0": [2, 3, 4], "X1": [2, 3, 4], "X2": [2, 3, 4], "XD": [2, 3, 4], "G1": [1], "G2": [2], "G3": [3],
      "X1n": [2], "C8": [2]}
MAIN = ("X0", "X1", "X2", "XD")
TARGETED = ("G1", "G2", "G3")
ORDER_MUTANTS = ("M1", "M2", "M3", "M5", "M6", "M7", "M8", "M15")
MIN_N = {"M1": 3, "M2": 3, "M3": 3, "M5": 3, "M6": 2, "M7": 2, "M8": 2, "M15": 2}
UNSUPPORTED_P6 = ("M4", "M4L", "M4R", "M9", "M10", "M11", "M12", "M13", "M14", "M16")


def on_path(m, N):
    return N >= MIN_N[m] and not (m == "M8" and N == 2)


# P5 recorded preflight status list (§2.4.3.6): not_infeasible at these five points, infeasible at the other 38.
P5_NOT_INFEASIBLE = frozenset({F(1, 2**19), F(3, 2**21), F(1, 2**20), F(3, 2**22), F(1, 2**21)})

# P7 evaluated sets
P7_SETS = {"empty": [], "minus_2^-40": [s for s in S_TEST if s != F(1, 2**40)],
           "plus_2^-41": S_TEST + [F(1, 2**41)],
           "replace_2^-19": [s for s in S_TEST if s != F(1, 2**19)] + [F(1, 2**19) + F(1, 2**60)]}

# T6.6 interface case data (§2.4.3.6)
K = [[F(3), F(0)], [F(4), F(0)]]
G_IFACE = 3 * 10**10

# C4 frozen inputs (§2.4.3.8)
C4_OPEN = [F(9, 512), F(5, 256)]
C4_CLOSED = [F(533, 512), F(1069, 1024)]
C4_EDGE = [F(1023, 1024), F(1025, 1024)]
C4_S = [F(1), F(1, 2**19), F(1, 2**40)]


def V4(j):
    return [[F((3 * a + 2 * j) % 9 - 4, 2) if a == b else
             F((-1) ** (a + b + j) * (1 + (min(a, b) + 2 * max(a, b) + j) % 3), 2) for b in range(4)] for a in range(4)]


def c4_lattices():
    out = []
    for g in C4_OPEN + C4_CLOSED + C4_EDGE:
        out.append((f"c4_ny1_{g}", lat([g], VS[:3]), [1, 2, 3]))
    for pair in ((C4_OPEN[0], C4_CLOSED[0]), (C4_OPEN[1], C4_CLOSED[1]), (C4_EDGE[0], C4_EDGE[1])):
        out.append((f"c4_ny2_{pair[0]}_{pair[1]}", lat(list(pair), VR[:3]), [1, 2, 3]))
    out.append(("c4_ny4_N6", lat(C4_OPEN + C4_CLOSED, [V4(j) for j in range(1, 7)]), [6]))
    return out
