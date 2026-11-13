"""T5 - exact-rational assembly (§2.4.2): (a) the binary64 import boundary, (b) L / R_W / R_orig rational assembly
against an independently written Fraction-only routine, float rejection, no float after import, the exactness check."""
import math
import pathlib
import sys
import types
from fractions import Fraction as F

import numpy as np

from . import inputs as I, indep
from .lattice import E_perp, _ROOT  # noqa: F401  (sys.path set for stationary)
from .result import Results
from stationary import Device, device_hamiltonian_blocks, scattering_matrix, transverse_basis  # noqa: E402

DX = 1.0 / 8


def _pow2(d):
    return d > 0 and d & (d - 1) == 0


def _float_tracer(prod_file, bad, seen=None):
    numeric = (float, complex, np.floating, np.complexfloating, np.ndarray)

    def tracer(frame, event, arg):
        if frame.f_code.co_filename != prod_file:
            return tracer
        if seen is not None:
            seen.append(frame.f_code.co_name)
        if event in ("line", "return"):
            for k, v in frame.f_locals.items():
                vals = v if isinstance(v, tuple) and len(v) <= 4 else (v,)
                for x in vals:
                    if isinstance(x, numeric) and not isinstance(x, bool):
                        bad.append((frame.f_code.co_name, k, type(x).__name__))
        return tracer
    return tracer


def _raises(fn, exc=TypeError):
    try:
        fn()
    except exc:
        return True
    except Exception:  # noqa: BLE001
        return False
    return False


def _is_binary64_format(t):
    """Independent test-side decision from the actual dtype: itemsize (per component) 8 bytes and a 52-bit mantissa."""
    f = np.finfo(t)
    comp = np.dtype(t).itemsize // (2 if np.issubdtype(t, np.complexfloating) else 1)
    return comp == 8 and f.nmant == 52 and f.nexp == 11


def _import_type_kats(R, prod):
    """Import-boundary type KATs (hardening of §2.4.1 / T5): accepted types convert exactly; every other precision
    format or object is rejected (never promoted, never narrowed)."""
    imp, impc = prod.import_binary64, prod.import_complex128
    third = F(6004799503160661, 2**54)  # binary64 nearest to 1/3
    # Python float and numpy.float64: accepted, exact
    R.add("T5 import: Python float accepted exactly", imp(1 / 3) == third, "CERT-VIOLATION")
    R.add("T5 import: numpy.float64 accepted exactly", imp(np.float64(1) / 3) == third, "CERT-VIOLATION")
    # lower precision: rejected (not promoted)
    for t in (np.float16, np.float32):
        R.add(f"T5 import: numpy.{t.__name__} rejected (lower precision, not promoted)",
              _raises(lambda: imp(t(0.5)), TypeError), "CERT-VIOLATION")
    # longdouble: platform dependent, decided from its actual format
    ld = np.longdouble(1) / 3
    if _is_binary64_format(np.longdouble):
        R.add(f"T5 import: numpy.longdouble is binary64 on this platform (itemsize {np.dtype(np.longdouble).itemsize}) "
              f"-> accepted exactly", imp(ld) == third, "CERT-VIOLATION")
    else:
        R.add(f"T5 import: numpy.longdouble is extended here (itemsize {np.dtype(np.longdouble).itemsize}, "
              f"{np.finfo(np.longdouble).nmant}-bit mantissa) -> rejected, not narrowed",
              _raises(lambda: imp(ld), TypeError), "CERT-VIOLATION")
    # Python complex and numpy.complex128: accepted, exact components
    want = (F(1, 2), -F(3602879701896397, 2**55))
    R.add("T5 import: Python complex accepted exactly", impc(0.5 - 0.1j) == want, "CERT-VIOLATION")
    R.add("T5 import: numpy.complex128 accepted exactly", impc(np.complex128(0.5 - 0.1j)) == want, "CERT-VIOLATION")
    R.add("T5 import: numpy.complex64 rejected (lower precision, not promoted)",
          _raises(lambda: impc(np.complex64(0.5 - 0.25j)), TypeError), "CERT-VIOLATION")
    cl = np.clongdouble(1) / 3 + 1j * (np.clongdouble(1) / 7)
    if _is_binary64_format(np.clongdouble):
        R.add(f"T5 import: numpy.clongdouble has binary64 components here (itemsize {np.dtype(np.clongdouble).itemsize})"
              f" -> accepted exactly", impc(cl) == (third, F(np.float64(1) / 7)), "CERT-VIOLATION")
    else:
        R.add(f"T5 import: numpy.clongdouble is extended here (itemsize {np.dtype(np.clongdouble).itemsize}, "
              f"{np.finfo(np.clongdouble).nmant}-bit mantissa) -> rejected, not narrowed",
              _raises(lambda: impc(cl), TypeError), "CERT-VIOLATION")
    # NaN / inf in either complex component
    for z, w in ((complex(float("nan"), 0), "NaN real part"), (complex(0, float("inf")), "+inf imaginary part"),
                 (np.complex128(complex(float("-inf"), 1)), "-inf real part (complex128)")):
        R.add(f"T5 import: complex with {w} rejected", _raises(lambda: impc(z), ValueError), "CERT-VIOLATION")

    class HasRealImag:  # an arbitrary object exposing .real and .imag
        real, imag = 0.5, 0.25

    class FloatSub(float):
        pass
    for obj, w, fn in ((HasRealImag(), "object with .real/.imag", impc), (HasRealImag(), "object with .real/.imag", imp),
                       (FloatSub(0.5), "float subclass", imp), (True, "bool", imp), (1, "int", imp), (F(1, 2), "Fraction", imp),
                       (1.0, "Python float passed to the complex importer", impc), (np.float64(1), "numpy.float64 to complex importer", impc),
                       (0.5 + 0j, "Python complex passed to the real importer", imp)):
        R.add(f"T5 import: {fn.__name__} rejects {w}", _raises(lambda: fn(obj), TypeError), "CERT-VIOLATION")
    # matrix import: float32 / complex64 / arbitrary entries rejected; float64 / complex128 accepted
    R.add("T5 import_matrix: float64 and complex128 arrays accepted exactly",
          prod.import_matrix(np.array([[0.5, 0.25]])) == (((F(1, 2), F(0)), (F(1, 4), F(0))),)
          and prod.import_matrix(np.array([[0.5 - 0.25j]])) == (((F(1, 2), F(-1, 4)),),), "CERT-VIOLATION")
    for arr, w in ((np.array([[0.5]], dtype=np.float32), "float32 array"), (np.array([[0.5]], dtype=np.complex64), "complex64 array"),
                   ([[HasRealImag()]], "arbitrary object entry")):
        R.add(f"T5 import_matrix rejects {w}", _raises(lambda: prod.import_matrix(arr), TypeError), "CERT-VIOLATION")


def run(prod):
    R = Results("T5 exact-rational assembly")
    # (a) import boundary
    imp = prod.import_binary64
    cases = [(0.1, F(3602879701896397, 2**55)), (-0.0, F(0)), (5e-324, F(1, 2**1074)),
             (1.7976931348623157e308, F(2**53 - 1) * F(2) ** 971), (math.pi, F(884279719003555, 2**48))]
    for x, exact in cases:
        q = imp(x)
        R.add(f"T5 import {x!r}: exact value, round trip, power-of-two denominator",
              q == exact and float(q) == x and _pow2(q.denominator) and isinstance(q, F), "CERT-VIOLATION", q)
    z = np.complex128(1.5 - 0.1j)
    qz = prod.import_complex128(z)
    R.add("T5 import complex128 (1.5 - 0.1j): exact parts and round trip",
          qz == (F(3, 2), -F(3602879701896397, 2**55)) and complex(float(qz[0]), float(qz[1])) == z, "CERT-VIOLATION", qz)
    for bad, w in ((float("nan"), "NaN"), (float("inf"), "+inf"), (float("-inf"), "-inf")):
        R.add(f"T5 import rejects {w}", _raises(lambda: imp(bad), ValueError), "CERT-VIOLATION")
    _import_type_kats(R, prod)
    # (b) rational assembly on three fixed configurations
    prod_file = str(pathlib.Path(prod.__file__).resolve())
    # tracer liveness controls: it must see production frames, and it must flag a planted float local
    def _planted(x):
        y = 0.5 * x
        return y

    planted = {"planted": types.FunctionType(_planted.__code__.replace(co_filename=prod_file, co_name="planted"), {})}
    bad0, seen0 = [], []
    old = sys.gettrace()
    sys.settrace(_float_tracer(prod_file, bad0, seen0))
    try:
        planted["planted"](F(1))
        prod.assemble_A1([(((F(1), F(0)),),)], "full")
    finally:
        sys.settrace(old)
    R.add("T5 tracer control: a planted float local in a production-file frame is flagged",
          any(b[0] == "planted" and b[2] == "float" for b in bad0), "KAT-FAIL", bad0[:3])
    R.add("T5 tracer control: production assembly frames are actually observed",
          "assemble_A1" in seen0 and "_mat" in seen0, "KAT-FAIL", sorted(set(seen0))[:6])
    for name, ny, N, cls, E, Jr, Ji in I.T5_CONFIGS:
        V1 = I.V1_blocks(ny, cls, Jr, Ji)
        Ep = E_perp(prod, ny)
        Epd = np.asarray(transverse_basis(1.0, ny).E_perp, dtype=float)
        V1d = np.array([[[complex(float(z[0]), float(z[1])) for z in r] for r in B] for B in V1])
        if cls == "R":
            V1d = V1d.real.copy()
        keys = ["full"] + [("prefix", j) for j in range(1, N)]
        A_L = {}
        for s in I.T5_S:
            sd = float(s)
            assert F(sd) == s
            Vs = sd * V1d
            Delta = (float(E) * np.eye(ny) - np.diag(Epd) - Vs) * (DX * DX)          # production W/Z input (binary64)
            dlead = (float(E) - Epd) * (DX * DX)
            dev = Device(xs=DX * np.arange(N), dx=DX, basis=transverse_basis(1.0, ny), V=Vs, V_lead_max=0.0)
            H, tB = device_hamiltonian_blocks(dev)
            A_orig = np.array([float(E) * np.eye(ny, dtype=complex) - H[j] for j in range(N)])
            sigma = np.asarray(scattering_matrix(float(E), dev, min_closed=0)["modes"]["sigma"])
            # import boundary: every binary64 input converted once
            Dq = [prod.import_matrix(Delta[j]) for j in range(N)]
            dq = [imp(float(x)) for x in dlead]
            Aq = [prod.import_matrix(A_orig[j]) for j in range(N)]
            sq = [prod.import_complex128(complex(x)) for x in sigma]
            for key in keys:
                kind, j = ("full", None) if key == "full" else ("prefix", key[1])
                n = N if key == "full" else key[1]
                bad, seen = [], []
                old = sys.gettrace()
                sys.settrace(_float_tracer(prod_file, bad, seen))
                try:
                    L = prod.assemble_L(E, Ep, V1, s, kind, j)
                    W = prod.assemble_R_W(Dq, dq, kind, j)
                    O = prod.assemble_R_orig(Aq, sq, kind, j)
                finally:
                    sys.settrace(old)
                tag = f"T5 {name} s={s} {key}"
                R.add(f"{tag}: no float/complex/numpy object in any production assembly frame after import "
                      f"({len(seen)} traced events)", not bad and {"assemble_L", "assemble_R_W", "assemble_R_orig"} <= set(seen),
                      "CERT-VIOLATION", bad[:3])
                R.add(f"{tag}: L == independent Fraction-only assembly",
                      L == indep.second_assembly_L(E, Ep, V1, s, n, key == "full"), "CERT-VIOLATION")
                R.add(f"{tag}: R_W == independent Fraction-only assembly",
                      W == indep.second_assembly_R_W(Dq, dq, n, key == "full"), "CERT-VIOLATION")
                R.add(f"{tag}: R_orig == independent Fraction-only assembly",
                      O == indep.second_assembly_R_orig(Aq, sq, n, key == "full"), "CERT-VIOLATION")
                A_L[(s, key)] = L
        for key in keys:  # exactness check: A_mid(2^-40) - A_mid(0) = -2^-40 blockdiag(V(1))
            n = N if key == "full" else key[1]
            D = tuple(tuple((x[0] - y[0], x[1] - y[1]) for x, y in zip(ra, rb))
                      for ra, rb in zip(A_L[(F(1, 2**40), key)], A_L[(F(0), key)]))
            want = [[(F(0), F(0))] * (n * ny) for _ in range(n * ny)]
            for b in range(n):
                for p in range(ny):
                    for r in range(ny):
                        want[b * ny + p][b * ny + r] = (-F(1, 2**40) * V1[b][p][r][0], -F(1, 2**40) * V1[b][p][r][1])
            R.add(f"T5 {name} {key}: A_mid(2^-40) - A_mid(0) = -2^-40 blockdiag(V(1)) exactly",
                  D == tuple(tuple(r) for r in want), "CERT-VIOLATION")
    # direct float input must be rejected by every rational assembly function (even exactly representable 0.5)
    name, ny, N, cls, E, Jr, Ji = I.T5_CONFIGS[1]
    V1 = I.V1_blocks(ny, cls, Jr, Ji)
    Ep = E_perp(prod, ny)
    Vf = [tuple(tuple((0.5, 0.0) for _ in r) for r in B) for B in V1]
    Dq0 = [tuple(tuple((F(0), F(0)) for _ in range(ny)) for _ in range(ny))] * N
    for desc, fn in (("assemble_L(E = 0.5)", lambda: prod.assemble_L(0.5, Ep, V1, F(0), "full")),
                     ("assemble_L(s = 0.5)", lambda: prod.assemble_L(E, Ep, V1, 0.5, "full")),
                     ("assemble_L(E_perp float)", lambda: prod.assemble_L(E, [float(e) for e in Ep], V1, F(0), "full")),
                     ("assemble_L(V float)", lambda: prod.assemble_L(E, Ep, Vf, F(0), "full")),
                     ("assemble_A1(V float)", lambda: prod.assemble_A1(Vf, "full")),
                     ("assemble_R_W(Delta float)", lambda: prod.assemble_R_W(Vf, [F(0)] * ny, "full")),
                     ("assemble_R_W(dlead float)", lambda: prod.assemble_R_W(Dq0, [0.5] * ny, "full")),
                     ("assemble_R_orig(A_orig float)", lambda: prod.assemble_R_orig(Vf, [(F(0), F(0))] * ny, "full")),
                     ("assemble_R_orig(sigma float)", lambda: prod.assemble_R_orig(Dq0, [(0.5, 0.0)] * ny, "full")),
                     ("assemble_A_mid(Delta float)", lambda: prod.assemble_A_mid(Vf, [F(1, 2)] * ny, "full")),
                     ("eps_A_L(E float)", lambda: prod.eps_A_L(0.5, Ep, N, "full")),
                     ("import_binary64(bool)", lambda: prod.import_binary64(True)),
                     ("import_binary64(int)", lambda: prod.import_binary64(3))):
        R.add(f"T5 type check: {desc} raises TypeError", _raises(fn, TypeError), "CERT-VIOLATION")
    return R
