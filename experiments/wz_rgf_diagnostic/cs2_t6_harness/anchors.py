"""Double-precision anchors (§2.4.3.5): reference exact endpoint blocks vs wz_rgf.wz_endpoint_blocks in binary64.
Tolerance: the §4.7 item 2 convention, ||X_double - X_exact||_F / ||G_hat||_2 <= 10 * cond_2(M) * eps.
DIAGNOSTIC ONLY: reported, never part of PASS; a disagreement only triggers an investigation note."""
import numpy as np

import wz_rgf
from . import inputs as I
from .common import ref_blocks, true_lams
from cs2_t6_reference import rlattice

EPS = np.finfo(float).eps


def _c(z):
    return complex(float(z.re), float(z.im))


def _np(X):
    return np.array([[_c(z) for z in r] for r in X], dtype=complex)


def run():
    rows = []
    for name in I.MAIN:
        lat = I.LATTICES[name]
        for N in I.NS[name]:
            for s in I.S6:
                b = ref_blocks(name, N, s)
                D = np.array([_np(Dj) for Dj in rlattice.delta_blocks(lat, N, s)])
                lm1 = np.array([_c(l) - 1.0 for l in true_lams(name)], dtype=complex)
                out = wz_rgf.wz_endpoint_blocks(D, lm1)
                Gh = _np(b["Ghat"])
                M = _np(b["A"]) / float(rlattice.T)
                cond = np.linalg.cond(M, 2)
                g2 = np.linalg.norm(Gh, 2)
                tol = 10 * cond * EPS
                worst = 0.0
                for k, rk in (("Ghat11", "G_11"), ("GhatN1", "G_N1"), ("Ghat1N", "G_1N"), ("GhatNN", "G_NN")):
                    worst = max(worst, np.linalg.norm(out[k] - _np(b[rk])) / g2)
                rows.append(dict(case=f"{name} N={N} s={s}", rel_err=float(worst), tol=float(tol), within=bool(worst <= tol)))
    return rows
