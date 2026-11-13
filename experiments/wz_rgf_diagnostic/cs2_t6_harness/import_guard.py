"""Stage 1: AST import rules of §2.4.3.5 (extended to the covering module and the T1-T5 package), checked before any
test runs. Dynamic-import bypasses (__import__, importlib.import_module, exec, eval) are forbidden in production and
reference code. The guard is self-tested on synthetic sources before it is applied."""
import ast
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
REF = ROOT / "cs2_t6_reference"
T1T5 = ROOT / "cs2_t1t5"
PROD = [ROOT / "coupling_scaling_v2_certify.py", ROOT / "coupling_scaling_v2_covering.py"]
HARN = HERE
REF_ALLOWED = {"fractions", "math", "functools", "__future__"}
REF_ALLOWED_RNUM = REF_ALLOWED | {"mpmath"}
PROD_ALLOWED = {"fractions", "math", "__future__", "mpmath", "numpy", "coupling_scaling_v2_certify"}
T1T5_FORBIDDEN = {"cs2_t6_harness", "wz_mutants", "coupled_blocks", "wz_rgf"}
DYNAMIC = {"__import__", "import_module", "exec", "eval"}


def _imports_src(src, name="<src>"):
    tree = ast.parse(src, filename=name)
    mods, dyn = [], []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods += [(a.name.split(".")[0], 0) for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            mods.append(((node.module or "").split(".")[0], node.level))
        elif isinstance(node, ast.Call):
            fn = node.func
            nm = fn.id if isinstance(fn, ast.Name) else fn.attr if isinstance(fn, ast.Attribute) else None
            if nm in DYNAMIC:
                dyn.append(nm)
    return mods, dyn


def _imports(path):
    return _imports_src(path.read_text(encoding="utf-8"), str(path))


def _ref_violations(src, allowed):
    mods, dyn = _imports_src(src)
    return [m for m, lvl in mods if lvl == 0 and m not in allowed] + dyn


def check(R):
    # self-test of the guard on synthetic sources
    R.add("guard control: forbidden import detected (import coupling_scaling_v2_certify)",
          _ref_violations("import coupling_scaling_v2_certify\n", REF_ALLOWED) == ["coupling_scaling_v2_certify"],
          "IMPORT-VIOLATION")
    R.add("guard control: forbidden from-import detected (from wz_mutants import x)",
          _ref_violations("from wz_mutants import x\n", REF_ALLOWED) == ["wz_mutants"], "IMPORT-VIOLATION")
    R.add("guard control: numpy detected in reference", _ref_violations("import numpy as np\n", REF_ALLOWED) == ["numpy"],
          "IMPORT-VIOLATION")
    R.add("guard control: dynamic __import__ / importlib.import_module / exec / eval detected",
          sorted(_ref_violations("__import__('os')\nimport importlib\nimportlib.import_module('x')\nexec('1')\neval('1')\n",
                                 REF_ALLOWED | {"importlib"})) == sorted(["__import__", "import_module", "exec", "eval"]),
          "IMPORT-VIOLATION")
    R.add("guard control: allowed reference imports pass", _ref_violations("from fractions import Fraction\nfrom . import rmat\n",
                                                                           REF_ALLOWED) == [], "IMPORT-VIOLATION")
    # reference
    for f in sorted(REF.glob("*.py")):
        allowed = REF_ALLOWED_RNUM if f.name == "rnum.py" else REF_ALLOWED
        mods, dyn = _imports(f)
        bad = [m for m, lvl in mods if lvl == 0 and m not in allowed] + dyn
        R.add(f"reference {f.name}: only {sorted(allowed)} (plus intra-package) imported; no dynamic import", not bad,
              "IMPORT-VIOLATION", f"forbidden: {bad}")
    mp_users = [f.name for f in sorted(REF.glob("*.py")) if any(m == "mpmath" for m, _ in _imports(f)[0])]
    R.add("only reference rnum.py imports mpmath", mp_users == ["rnum.py"], "IMPORT-VIOLATION", mp_users)
    # production
    for f in PROD:
        mods, dyn = _imports(f)
        bad = [m for m, lvl in mods if m not in PROD_ALLOWED or m.startswith("cs2_")] + dyn
        R.add(f"production {f.name}: no reference/harness/T1-T5/repository-solver import; no dynamic import", not bad,
              "IMPORT-VIOLATION", f"forbidden: {bad}")
    # T1-T5 package (system under test = production; references written there or clean-room exact algebra)
    for f in sorted(T1T5.glob("*.py")):
        mods, dyn = _imports(f)
        bad = [m for m, lvl in mods if lvl == 0 and m in T1T5_FORBIDDEN] + dyn
        R.add(f"T1-T5 {f.name}: no harness / W/Z solver import; no dynamic import", not bad, "IMPORT-VIOLATION", f"{bad}")
    # harness: repository solvers only in anchors.py
    for f in sorted(HARN.glob("*.py")):
        if f.name in ("anchors.py", "import_guard.py"):
            continue
        mods, _ = _imports(f)
        bad = [m for m, lvl in mods if lvl == 0 and m in {"wz_mutants", "wz_rgf", "coupled_blocks", "stationary"}]
        R.add(f"harness {f.name}: repository solvers only in anchors.py", not bad, "IMPORT-VIOLATION", f"{bad}")
    return R
