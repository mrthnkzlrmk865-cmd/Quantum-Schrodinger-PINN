"""L-lattice physical inputs: E (exact dyadic), E_perp imported from stationary.transverse_basis (binary64 -> exact)."""
from fractions import Fraction as F
from functools import lru_cache

import pathlib
import sys

_ROOT = str(pathlib.Path(__file__).resolve().parents[3])   # repository root (stationary.py), as coupled_blocks.py does
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from stationary import transverse_basis  # noqa: E402  (read-only use: E_perp source values)
from . import inputs as I


@lru_cache(maxsize=None)
def E_perp(prod, ny):
    return tuple(prod.import_binary64(float(e)) for e in transverse_basis(1.0, ny).E_perp)


def matrices(prod, cfg, s):
    """Production A_mid(s), A1 and eps_A for the full matrix and every prefix of an L configuration."""
    name, ny, N, cls, E, Jr, Ji = cfg
    V1 = I.V1_blocks(ny, cls, Jr, Ji)
    Ep = E_perp(prod, ny)
    out = {}
    for key in ["full"] + [("prefix", j) for j in range(1, N)]:
        kind, j = ("full", None) if key == "full" else ("prefix", key[1])
        out[key] = dict(A=prod.assemble_L(E, Ep, V1, s, kind, j), A1=prod.assemble_A1(V1, kind, j),
                        eps=prod.eps_A_L(E, Ep, N, kind), n=(N if key == "full" else key[1]), full=(key == "full"))
    return out, V1, Ep
