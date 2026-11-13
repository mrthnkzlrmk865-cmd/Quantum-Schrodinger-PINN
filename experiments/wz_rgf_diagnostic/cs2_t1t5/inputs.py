"""Literal T1-T5 test-file inputs (§2.4.2 delegates the fixed integers to the test file). Fixed before any result."""
from fractions import Fraction as F

T = F(32)
DX2 = F(1, 64)
V0 = 8
T5_S = [F(0), F(1), F(3, 4), F(1, 2**40), F(3, 2**42)]

# Jr_j / Ji_j integer arrays in [-8, 8] (explicit), and explicit dyadic energies.
JR_NY2 = [[[3, -5], [2, 7]], [[-6, 1], [4, -2]], [[8, -3], [-7, 5]], [[-1, 6], [3, -4]], [[2, 2], [-8, 1]], [[-5, 4], [6, -3]]]
JI_NY2 = [[[0] * 2 for _ in range(2)] for _ in range(6)]
JR_NY3 = [[[1, -4, 6], [7, -2, 3], [-5, 8, 0]], [[-3, 5, -1], [2, 4, -6], [8, -7, 1]]]
JI_NY3 = [[[2, 0, -3], [5, -1, 4], [-6, 7, 2]], [[-4, 3, 1], [0, 6, -5], [3, -2, -8]]]


def _jr4(j, a, b):
    return ((3 * a + 5 * b + 7 * j) % 17) - 8          # explicit deterministic integers in [-8, 8]


def _ji4(j, a, b):
    return ((5 * a + 2 * b + 3 * j) % 15) - 7


JR_NY4 = [[[_jr4(j, a, b) for b in range(4)] for a in range(4)] for j in range(6)]
JI_NY4 = [[[_ji4(j, a, b) for b in range(4)] for a in range(4)] for j in range(6)]

# configurations: (name, N_y, N, class, E, Jr, Ji)
T5_CONFIGS = [("t5_ny2_N1", 2, 1, "R", F(640, 64), JR_NY2[:1], JI_NY2[:1]),
              ("t5_ny3_N2", 3, 2, "H", F(1600, 64), JR_NY3, JI_NY3),
              ("t5_ny4_N6", 4, 6, "H", F(3200, 64), JR_NY4, JI_NY4)]
T1F_CONFIGS = [("t1f_ny2_N1", 2, 1, "R", F(640, 64), JR_NY2[:1], JI_NY2[:1]),
               ("t1f_ny2_N2", 2, 2, "R", F(640, 64), JR_NY2[:2], JI_NY2[:2]),
               ("t1f_ny2_N6", 2, 6, "R", F(640, 64), JR_NY2, JI_NY2),
               ("t1f_ny4_N6", 4, 6, "H", F(3200, 64), JR_NY4, JI_NY4)]
# T4: N = 2, N_y = 2 L-lattice with V_1(1) = 0 and V_2(1) != 0 (fixed integers)
T4_CONFIG = ("t4_ny2_N2", 2, 2, "R", F(640, 64), [[[0, 0], [0, 0]], JR_NY2[1]], [[[0, 0], [0, 0]], [[0, 0], [0, 0]]])

# T3: exact rational deltas, and the fixed N_y = 4 configuration (explicit integers / 64)
T3_DELTAS = [F(1, 50), F(-1, 50), F(99, 50), F(1), F(1, 2), F(-1, 2), F(-1, 4)]
T3_NY4_DELTAS = [F(75, 64), F(40, 64), F(-13, 64), F(-90, 64)]


def V1_blocks(ny, cls, Jr, Ji):
    """V_j(1) = v0 (B_j + B_j^H)/2 with B_j = (Jr_j + i Ji_j [class H]) / 8  (§2.3), as exact (re, im) pairs."""
    out = []
    for jr, ji in zip(Jr, Ji):
        blk = []
        for a in range(ny):
            row = []
            for b in range(ny):
                re = F(V0, 8) * F(jr[a][b] + jr[b][a], 2)
                im = F(V0, 8) * F(ji[a][b] - ji[b][a], 2) if cls == "H" else F(0)
                if cls == "D" and a != b:
                    re = im = F(0)
                row.append((re, im))
            blk.append(tuple(row))
        out.append(tuple(blk))
    return out
