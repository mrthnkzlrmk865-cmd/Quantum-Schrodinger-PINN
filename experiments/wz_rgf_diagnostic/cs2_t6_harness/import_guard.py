"""Stage 1: AST import rules of §2.4.3.5, checked before any test runs."""
import ast
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
REF = ROOT / "cs2_t6_reference"
PROD = ROOT / "coupling_scaling_v2_certify.py"
HARN = HERE
REPO_MODULES = {"coupling_scaling_v2_certify", "wz_mutants", "coupled_blocks", "wz_rgf", "stationary", "wz_variants",
                "cs2_t6_harness", "smatrix_spec", "sparse_reference", "config", "transport"}
REF_ALLOWED = {"fractions", "math", "functools", "__future__"}
REF_ALLOWED_RNUM = REF_ALLOWED | {"mpmath"}


def _imports(path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out += [(a.name.split(".")[0], 0) for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            out.append(((node.module or "").split(".")[0], node.level))
        elif isinstance(node, ast.Call) and getattr(node.func, "id", "") in ("__import__",):
            out.append(("<dynamic __import__>", 0))
        elif isinstance(node, ast.Attribute) and node.attr == "import_module":
            out.append(("<importlib.import_module>", 0))
    return out


def check(R):
    for f in sorted(REF.glob("*.py")):
        allowed = REF_ALLOWED_RNUM if f.name == "rnum.py" else REF_ALLOWED
        bad = [m for m, lvl in _imports(f) if lvl == 0 and m not in allowed]
        R.add(f"reference {f.name}: only {sorted(allowed)} (plus intra-package) imported", not bad, "IMPORT-VIOLATION",
              f"forbidden: {bad}")
    bad = [m for m, lvl in _imports(PROD) if m.startswith("cs2_t6") or m.startswith("<")]
    R.add("production imports neither cs2_t6_reference nor cs2_t6_harness (no dynamic imports)", not bad,
          "IMPORT-VIOLATION", f"forbidden: {bad}")
    allowed_anchor = {"anchors.py"}
    for f in sorted(HARN.glob("*.py")):
        if f.name in allowed_anchor or f.name == "import_guard.py":
            continue
        bad = [m for m, lvl in _imports(f) if lvl == 0 and m in {"wz_mutants", "wz_rgf", "coupled_blocks", "stationary"}]
        R.add(f"harness {f.name}: repository solvers only in anchors.py", not bad, "IMPORT-VIOLATION", f"{bad}")
    return R
