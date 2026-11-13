"""Stage 6: mutation suite MU1-MU27 (§2.4.3.10). Each mutation is an exact-once source patch applied to a COPY of the
production module (never written to disk). Detection = at least one CERT-VIOLATION / PATH-FAIL / INTERFACE-FAIL finding,
or an exception raised by the mutant, in the listed detecting test(s). PRECONDITION-UNMET never counts as detection."""
import pathlib
import types
from fractions import Fraction as F

from . import inputs as I
from .common import Results, all_cert_cases
from . import t6_1, t6_2, t6_3, t6_4, t6_5, t6_6, t6_7, t6_8

PROD_PATH = pathlib.Path(__file__).resolve().parent.parent / "coupling_scaling_v2_certify.py"
DETECTING = {"CERT-VIOLATION", "PATH-FAIL", "INTERFACE-FAIL", "MUTANT-RAISED"}
TARGETED_CASES = [c for c in all_cert_cases() if c[0] in I.TARGETED]
X12_CASES = [c for c in all_cert_cases() if c[0] in ("X1", "X2")]


def _swap(line, var):
    indent = line[:len(line) - len(line.lstrip())]
    return line + f"\n{indent}{var}[0], {var}[1] = {var}[1], {var}[0]"


L_M13 = '        longest = ["G_NN", "G1_11"] + [("P", j) for j in range(2, N)]  # M1/M3 longest'
L_M2 = '        longest = [("P", j) for j in range(N - 1, 1, -1)] + ["G1_11", "G_NN"]  # M2 longest'
L_M5 = '            first = [("chi", j - 1), ("P", j), ("psi", j - 1)]  # M5 longest'
L_M6 = '        longest = [("chi", N - 1), "G_NN"]  # M6/M15 longest'
L_M7 = '        longest = ["G_NN", ("psi", N - 1)]  # M7 longest'
L_M8 = '        longest = [("chi", N - 1), "G_NN", ("psi", N - 1)]  # M8 longest'
GLO = "    G_lo = max(T_HOP / h_W, T_HOP * _lb(maxabs2) - eta_full)"
TH = "THRESHOLD = Fraction(1, 10**10)"


def t64(cases):
    return lambda p: t6_4.run(p, cases=cases)


def t68(*paths):
    return lambda p: t6_8.run(p, only=set(paths))


MUTATIONS = [
    ("MU1", "eps_A full N>=2: sum instead of 2*sum",
     [("return root(2 * s_te2)                        # full N >= 2", "return root(s_te2)  # MU1")], [t64(TARGETED_CASES)]),
    ("MU2", "eps_A N=1: (te)^2 instead of (2te)^2",
     [("return root(s_2te2)                       # full N = 1", "return root(s_te2)  # MU2")], [t64(TARGETED_CASES)]),
    ("MU3", "eps_A N=1 uses the N>=2 formula",
     [("return root(s_2te2)                       # full N = 1", "return root(2 * s_te2)  # MU3")], [t64(TARGETED_CASES)]),
    ("MU4", "eps_A with lb in place of ub", [("root = _ub  # eps_A rounding", "root = _lb  # MU4")], [t64(TARGETED_CASES)]),
    ("MU5", "eps_A = 0 for prefixes",
     [("return root(s_te2)                        # prefix", "return Fraction(0)  # MU5")], [t64(TARGETED_CASES)]),
    ("MU6", "closed-branch midpoint gamma + mu",
     [('lam_mid = (g - mu, Fraction(0)) if branch == "closed"', 'lam_mid = (g + mu, Fraction(0)) if branch == "closed"')],
     [t64(X12_CASES)]),
    ("MU7.M1/M3", "swap first two factors of the longest monomial (M1/M3)", [(L_M13, _swap(L_M13, "longest"))], [t6_1.run]),
    ("MU7.M2", "swap first two factors of the longest monomial (M2)", [(L_M2, _swap(L_M2, "longest"))], [t6_1.run]),
    ("MU7.M5", "swap first two factors of each first M5 monomial", [(L_M5, _swap(L_M5, "first"))], [t6_1.run]),
    ("MU7.M6/M15", "swap first two factors of the longest monomial (M6/M15)", [(L_M6, _swap(L_M6, "longest"))], [t6_1.run]),
    ("MU7.M7", "swap first two factors of the longest monomial (M7)", [(L_M7, _swap(L_M7, "longest"))], [t6_1.run]),
    ("MU7.M8", "swap first two factors of the longest monomial (M8)", [(L_M8, _swap(L_M8, "longest"))], [t6_1.run]),
    ("MU8.M1/M3", "flip the single-block sign (M1/M3)",
     [('single = (-1, ["G_N1"])  # M1/M3 single', 'single = (+1, ["G_N1"])  # MU8')], [t6_1.run]),
    ("MU8.M2", "flip the single-block sign (M2)",
     [('single = (-1, ["G_1N"])  # M2 single', 'single = (+1, ["G_1N"])  # MU8')], [t6_1.run]),
    ("MU8.M6/M15", "flip the single-block sign (M6/M15)",
     [('single = (-1, ["G_N1"])  # M6/M15 single', 'single = (+1, ["G_N1"])  # MU8')], [t6_1.run]),
    ("MU8.M7", "flip the single-block sign (M7)",
     [('single = (-1, ["G_1N"])  # M7 single', 'single = (+1, ["G_1N"])  # MU8')], [t6_1.run]),
    ("MU9", "M5 sum over j = 2..N-2", [("for j in range(2, N):  # M5 terms", "for j in range(2, N - 1):  # MU9")], [t6_1.run]),
    ("MU10", "E_m := 0", [("E_m += abs(c) * (prod_up - prod)", "E_m += 0  # MU10")], [t6_2.run]),
    ("MU11", "U_m without the ub(N_y)*E_m term",
     [("U_m = _ub(_fro2(F)) + _ub(Fraction(N_y)) * E_m", "U_m = _ub(_fro2(F))  # MU11")], [t6_5.run]),
    ("MU12", "a := ub(f_Y) (no 1/(1 - omega))", [("a = _ub(f_Y) / (1 - omega)", "a = _ub(f_Y)  # MU12")], [t6_3.run]),
    ("MU13", "G_lo with +eta", [(GLO, GLO.replace("- eta_full)", "+ eta_full)  # MU13"))], [t6_6.run]),
    ("MU14", "no omega < 1 check", [("certified = omega < 1", "certified = True  # MU14")], [t68("P1", "P2", "P4")]),
    ("MU15", "THRESHOLD = Fraction(1e-10)", [(TH, "THRESHOLD = Fraction(1e-10)  # MU15")], [t6_7.run]),
    ("MU16", "THRESHOLD = 1e-10 (float)", [(TH, "THRESHOLD = 1e-10  # MU16")], [t6_7.run]),
    ("MU17", "U <= threshold*G_lo", [("return U < threshold * G_lo", "return U <= threshold * G_lo  # MU17")], [t6_7.run]),
    ("MU18", "THRESHOLD = Fraction(1, 10**9)", [(TH, "THRESHOLD = Fraction(1, 10**9)  # MU18")], [t6_7.run, t68("P5")]),
    ("MU19", "S_TEST with k = 19..39 only (42 points)",
     [("S_TEST = frozenset(s for s in _GRID if s <= Fraction(1, 2**19))",
       "S_TEST = frozenset(s for s in _GRID if s <= Fraction(1, 2**19) and s != Fraction(1, 2**40))  # MU19")],
     [t68("edge", "P7")]),
    ("MU20", "conclude: any point instead of every point",
     [('if all(r["point"] == "infeasible" for r in rs):', 'if any(r["point"] == "infeasible" for r in rs):  # MU20')],
     [t68("P5")]),
    ("MU21", "conclude skips the evaluated-set check",
     [("if frozenset(per_point) != S_TEST:", "if False:  # MU21")], [t68("P7")]),
    ("MU22", "unsupported mutant evaluated with the M6 polynomial",
     [('        return _inconclusive("unsupported_mutant")', '        mutant = "M6"  # MU22')], [t68("P6")]),
    ("MU23", "missing prefix replaced by zero matrices",
     [('            return _inconclusive("missing_reference")',
       '            refs = {**refs, key: (_zeros(key[1] * N_y), _zeros(key[1] * N_y), Fraction(0))}  # MU23')],
     [t68("P3")]),
    ("MU24", "point_infeasible without the G_lo < 0 rejection",
     [("if G_lo < 0:  # generic lower-bound check", "if False:  # MU24")], [t6_6.run]),
    ("MU25", "point_infeasible without the U < 0 rejection",
     [("if U < 0:  # generic upper-bound check", "if False:  # MU25")], [t6_6.run]),
    ("MU26", "bool accepted as an exact operand",
     [("    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):\n        raise TypeError(f\"{what}: exact int or",
       "    if not isinstance(x, (int, Fraction)):  # MU26\n        raise TypeError(f\"{what}: exact int or")], [t6_6.run]),
    ("MU27", "G_lo without the t/h_W branch", [(GLO, "    G_lo = T_HOP * _lb(maxabs2) - eta_full  # MU27")], [t6_6.run]),
]

# Supplementary (not part of MU1-MU27): exercises the N4 route directly, because the frozen MU18 constant change is also
# rejected by point_infeasible's literal threshold check (which raises) before the N4 mechanism can act.
SUPPLEMENTARY = [
    ("MU18-N4", "comparison uses 10*threshold (threshold misread by a decade, check bypassed)",
     [("return U < threshold * G_lo", "return U < 10 * threshold * G_lo  # MU18-N4")], [t68("P5")]),
]


def build_mutant(mid, patches):
    src = PROD_PATH.read_text(encoding="utf-8")
    for old, new in patches:
        n = src.count(old)
        if n != 1:
            raise RuntimeError(f"{mid}: anchor occurs {n} times (must be exactly once): {old[:70]!r}")
        src = src.replace(old, new)
    mod = types.ModuleType(f"coupling_scaling_v2_certify__{mid}")
    mod.__file__ = f"<mutant {mid}>"
    exec(compile(src, f"<mutant {mid}>", "exec"), mod.__dict__)
    return mod


def run_one(mid, desc, patches, detectors):
    mod = build_mutant(mid, patches)
    findings, mechanism = [], []
    for det in detectors:
        try:
            R = det(mod)
            hits = [f for f in R.failures if f["cls"] in DETECTING]
            findings += hits
            mechanism += sorted({f["cls"] for f in hits})
        except Exception as e:  # noqa: BLE001
            findings.append(dict(name=f"{mid} raised", ok=False, cls="MUTANT-RAISED", detail=f"{type(e).__name__}: {e}"))
            mechanism.append("MUTANT-RAISED")
    return dict(id=mid, description=desc, detected=bool(findings), n_findings=len(findings),
                mechanism=sorted(set(mechanism)), first=[f["name"] + " :: " + f["detail"][:120] for f in findings[:3]])


def run(include_supplementary=True):
    R = Results("6 mutation suite")
    rows = []
    for mid, desc, patches, dets in MUTATIONS + (SUPPLEMENTARY if include_supplementary else []):
        row = run_one(mid, desc, patches, dets)
        rows.append(row)
        R.add(f"{mid} detected ({desc})", row["detected"], "MUTATION-SURVIVED", "; ".join(row["mechanism"]))
    return R, rows
