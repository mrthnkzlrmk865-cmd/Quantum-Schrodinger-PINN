"""§5.6 polynomial forms F_m as monomial lists, their designated blocks, and evaluation on block dictionaries.

Keys: "G_NN", "G_N1", "G_1N", "G1_11", ("P", j), ("chi", j), ("psi", j). Used only after being proved equal to the exact
recursion difference (rmutants - rblocks) on every T6.1 input.
"""
from .qexact import RC
from . import rmat

SUPPORTED = ("M1", "M2", "M3", "M5", "M6", "M7", "M8", "M15")
MIN_N = {"M1": 3, "M2": 3, "M3": 3, "M5": 3, "M6": 2, "M7": 2, "M8": 2, "M15": 2}
DESIGNATED = {"M1": "G_N1", "M3": "G_N1", "M6": "G_N1", "M15": "G_N1", "M2": "G_1N", "M7": "G_1N",
              "M5": "G_11", "M8": "G_11"}


def on_path(m, N):
    return m in SUPPORTED and N >= MIN_N[m] and not (m == "M8" and N == 2)


def monomials(m, N):
    if m in ("M1", "M3"):
        return [(-1, ["G_NN", "G1_11"] + [("P", j) for j in range(2, N)]), (-1, ["G_N1"])]
    if m == "M2":
        return [(-1, [("P", j) for j in range(N - 1, 1, -1)] + ["G1_11", "G_NN"]), (-1, ["G_1N"])]
    if m == "M5":
        out = []
        for j in range(2, N):
            out.append((-1, [("chi", j - 1), ("P", j), ("psi", j - 1)]))
            out.append((+1, [("psi", j - 1), ("P", j), ("chi", j - 1)]))
        return out
    if m in ("M6", "M15"):
        return [(-1, [("chi", N - 1), "G_NN"]), (-1, ["G_N1"])]
    if m == "M7":
        return [(-1, ["G_NN", ("psi", N - 1)]), (-1, ["G_1N"])]
    if m == "M8":
        return [(+1, [("chi", N - 1), "G_NN", ("psi", N - 1)]), (-1, [("psi", N - 1), "G_NN", ("chi", N - 1)])]
    raise ValueError(f"unsupported mutant {m}")


def evaluate(monos, blocks):
    ny = len(blocks["G_NN"])
    F = rmat.zeros(ny)
    for c, keys in monos:
        X = blocks[keys[0]]
        for k in keys[1:]:
            X = rmat.matmul(X, blocks[k])
        F = rmat.add(F, rmat.scale(RC(c), X))
    return F


def F_ref(m, N, blocks):
    return evaluate(monomials(m, N), blocks)
