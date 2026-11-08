r"""
Fail-closed gate registry and overall-status aggregation for validation / physical campaigns.

Design (docs/PHYSICAL_CAMPAIGN_PLAN.md section 8):
  * The registry of MANDATORY gates is built from the predeclared protocol (expected row keys, thresholds), never from whatever rows a
    run happened to produce.
  * Every gate has exactly one state in STATES = (PASS, FAIL, NOT RUN, NOT EVALUABLE). Overall PASS iff every mandatory gate is
    explicitly PASS. A missing gate, an unknown state, an exception, a missing table, a missing expected row, a duplicate row, an
    unknown row status or a stage that did not complete all prevent PASS.
  * Statuses inside nested upstream summaries are flattened recursively, so a failure cannot hide in a sub-dictionary.
  * Rows are only read, never filtered out: failed and NOT EVALUABLE rows stay in the tables and are counted.

Two gate-specification builders are provided:
  b3y_gate_spec(protocol)                 from docs/rgf_b3y_protocol.json (read-only use: regression against existing B3Y results)
  physical_gate_spec(draft, n_seeds)      from docs/physical_campaign_gates_draft.json (DRAFT; not frozen)

This module does not import, and does not depend on, validate_rgf_b3y.py.
Command line (read-only):  python campaign_status.py --b3y results/rgf_b3y_validation
"""

from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass, field
from typing import Callable, Optional

PASS, FAIL, NOT_RUN, NOT_EVALUABLE = "PASS", "FAIL", "NOT RUN", "NOT EVALUABLE"
STATES = (PASS, FAIL, NOT_RUN, NOT_EVALUABLE)
OVERALL_INCOMPLETE = "INCOMPLETE"
HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------
# value normalisation
# ----------------------------------------------------------------------
def norm(v):
    """Canonical comparison value: numbers compare as floats (so '25', '25.0' and 25 are equal), everything else as stripped text."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, float)):
        return float(v) + 0.0
    s = str(v).strip()
    try:
        return float(s) + 0.0
    except ValueError:
        return s


def num(v):
    """Strict float parse; None if not a finite-or-infinite number (empty strings and text are not numbers)."""
    if isinstance(v, bool) or v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def is_true(v):
    return v is True or (isinstance(v, str) and v.strip() == "True")


# ----------------------------------------------------------------------
# checks and gate specifications
# ----------------------------------------------------------------------
@dataclass
class Check:
    column: str
    op: str                  # 'le' (value <= limit), 'eq' (numeric equality), 'true', 'finite', 'nonempty', 'eqcol' (== another column)
    limit: object = None

    def evaluate(self, row) -> Optional[str]:
        """None if the check passes, otherwise a reason. Missing columns and unparsable values fail (fail-closed)."""
        if self.column not in row:
            return f"column {self.column!r} missing"
        v = row[self.column]
        if self.op == "true":
            return None if is_true(v) else f"{self.column}={v!r} is not True"
        if self.op == "nonempty":
            return None if str(v).strip() else f"{self.column} empty"
        x = num(v)
        if x is None:
            return f"{self.column}={v!r} is not a number"
        if self.op == "finite":
            return None if math.isfinite(x) else f"{self.column}={x} not finite"
        if self.op == "le":
            return None if (math.isfinite(x) and x <= self.limit) else f"{self.column}={x!r} > {self.limit!r}"
        if self.op == "eq":
            return None if x == self.limit else f"{self.column}={x!r} != {self.limit!r}"
        if self.op == "eqcol":
            y = num(row.get(self.limit))
            return None if (y is not None and x == y) else f"{self.column}={x!r} != {self.limit}={row.get(self.limit)!r}"
        return f"unknown check operator {self.op!r}"


@dataclass
class GateSpec:
    gate_id: str
    description: str
    table: str                                         # file name inside the output directory
    key_columns: tuple
    expected_keys: Optional[list] = None               # protocol-derived row keys (tuples); every key must be present exactly once
    expected_count: Optional[int] = None               # alternative: number of rows (keys chosen at run time, e.g. a random subset)
    row_filter: dict = field(default_factory=dict)     # column -> value; selects the rows belonging to this gate
    checks: list = field(default_factory=list)
    status_column: Optional[str] = None
    status_ok: tuple = ()
    status_not_evaluable: tuple = ()
    allow_duplicates: bool = False
    row_evaluator: Optional[Callable] = None           # row -> None | reason (e.g. the TOST equivalence + power check)
    extra_rows_allowed: bool = True                    # extra rows are evaluated like expected rows, never ignored


def _row_key(row, cols):
    return tuple(norm(row.get(c)) for c in cols)


def read_table(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def evaluate_gate(spec: GateSpec, out_dir: str, cache: dict) -> dict:
    """Evaluate one mandatory gate. Never raises: any exception is a FAIL of that gate."""
    rec = dict(gate_id=spec.gate_id, description=spec.description, table=spec.table, state=FAIL, reason="", n_expected=None,
               n_rows=0, n_missing=0, n_duplicate_keys=0, n_failed_rows=0, n_not_evaluable_rows=0, n_extra_rows=0, first_failures="")
    try:
        path = os.path.join(out_dir, spec.table)
        if spec.table not in cache:
            cache[spec.table] = read_table(path) if os.path.exists(path) else None
        table = cache[spec.table]
        if table is None:
            rec.update(state=NOT_RUN, reason=f"table {spec.table} missing")
            return rec
        flt = {c: norm(v) for c, v in spec.row_filter.items()}
        rows = [r for r in table if all(norm(r.get(c)) == v for c, v in flt.items())]
        rec["n_rows"] = len(rows)
        reasons = []
        if spec.expected_keys is not None:
            exp = [tuple(norm(x) for x in k) for k in spec.expected_keys]
            rec["n_expected"] = len(exp)
            if len(set(exp)) != len(exp):
                reasons.append("protocol-derived expected keys are not unique")
            seen = {}
            for r in rows:
                seen.setdefault(_row_key(r, spec.key_columns), []).append(r)
            exp_set = set(exp)
            missing = [k for k in exp if k not in seen]
            dup = [k for k, v in seen.items() if len(v) > 1]
            extra = [k for k in seen if k not in exp_set]
            rec.update(n_missing=len(missing), n_duplicate_keys=len(dup), n_extra_rows=len(extra))
            if missing:
                reasons.append(f"{len(missing)} expected rows missing (e.g. {missing[0]})")
            if dup and not spec.allow_duplicates:
                reasons.append(f"{len(dup)} duplicated keys (e.g. {dup[0]})")
            if extra and not spec.extra_rows_allowed:
                reasons.append(f"{len(extra)} unexpected rows")
        if spec.expected_count is not None:
            rec["n_expected"] = spec.expected_count
            if len(rows) != spec.expected_count:
                reasons.append(f"{len(rows)} rows, protocol expects {spec.expected_count}")
        if rec["n_expected"] is None:
            reasons.append("gate has neither expected keys nor an expected count (protocol derivation missing)")
        failed, ne = [], 0
        for r in rows:
            if spec.status_column is not None:
                st = str(r.get(spec.status_column, "")).strip()
                if st in spec.status_not_evaluable:
                    ne += 1
                    continue
                if st not in spec.status_ok:
                    failed.append(f"{_row_key(r, spec.key_columns)}: status {st!r}")
                    continue
            for ch in spec.checks:
                why = ch.evaluate(r)
                if why:
                    failed.append(f"{_row_key(r, spec.key_columns)}: {why}")
                    break
            else:
                if spec.row_evaluator is not None:
                    why = spec.row_evaluator(r)
                    if why:
                        failed.append(f"{_row_key(r, spec.key_columns)}: {why}")
        rec["n_failed_rows"], rec["n_not_evaluable_rows"] = len(failed), ne
        if failed:
            reasons.append(f"{len(failed)} rows fail")
            rec["first_failures"] = " | ".join(failed[:3])
        if reasons:
            rec.update(state=FAIL, reason="; ".join(reasons))
        elif ne:
            rec.update(state=NOT_EVALUABLE, reason=f"{ne} NOT EVALUABLE rows")
        else:
            rec.update(state=PASS, reason="all expected rows present and passing")
    except Exception as e:                                                    # noqa: BLE001  (fail closed)
        rec.update(state=FAIL, reason=f"exception during evaluation: {type(e).__name__}: {e}")
    return rec


# ----------------------------------------------------------------------
# nested status flattening (upstream summaries)
# ----------------------------------------------------------------------
STATUS_KEYS = ("status", "state", "passed")


def flatten_statuses(obj, path=""):
    """All (path, value) pairs whose key is a status-like key, at any nesting depth (dicts and lists)."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{path}/{k}"
            if k in STATUS_KEYS and not isinstance(v, (dict, list)):
                out.append((p, v))
            else:
                out += flatten_statuses(v, p)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += flatten_statuses(v, f"{path}[{i}]")
    return out


def nested_status_gate(gate_id, summary_path, subtree, informational=("CLASSIFICATION ONLY",), require_null=("stopped_at",)) -> dict:
    """Gate on an upstream summary: every nested status must be PASS (or an explicitly informational label), every nested 'passed' True,
    and the early-termination fields must be null."""
    rec = dict(gate_id=gate_id, description=f"nested statuses of {os.path.basename(summary_path)}:{subtree}", table=os.path.basename(summary_path),
               state=FAIL, reason="", n_expected=None, n_rows=0, n_missing=0, n_duplicate_keys=0, n_failed_rows=0, n_not_evaluable_rows=0,
               n_extra_rows=0, first_failures="")
    try:
        if not os.path.exists(summary_path):
            rec.update(state=NOT_RUN, reason="summary missing")
            return rec
        with open(summary_path) as f:
            s = json.load(f)
        bad = []
        for k in require_null:
            if k not in s:
                bad.append(f"{k} missing")
            elif s[k] is not None:
                bad.append(f"{k} = {s[k]!r} (stopped early)")
        if subtree not in s:
            bad.append(f"{subtree} missing")
            flat = []
        else:
            flat = flatten_statuses(s[subtree], subtree)
        rec["n_rows"] = len(flat)
        if not flat and subtree in s:
            bad.append("no statuses found")
        for p, v in flat:
            if p.endswith("/passed"):
                if v is not True:
                    bad.append(f"{p} = {v!r}")
            elif v == PASS or v in informational:
                continue
            elif v == NOT_EVALUABLE:
                bad.append(f"{p} = NOT EVALUABLE")
            else:
                bad.append(f"{p} = {v!r}")
        rec["n_failed_rows"] = len(bad)
        if bad:
            ne_only = all("NOT EVALUABLE" in b for b in bad)
            rec.update(state=NOT_EVALUABLE if ne_only else FAIL, reason=f"{len(bad)} nested non-PASS statuses", first_failures=" | ".join(bad[:3]))
        else:
            rec.update(state=PASS, reason=f"{len(flat)} nested statuses PASS/informational")
    except Exception as e:                                                    # noqa: BLE001
        rec.update(state=FAIL, reason=f"exception: {type(e).__name__}: {e}")
    return rec


# ----------------------------------------------------------------------
# overall status
# ----------------------------------------------------------------------
def combine(registry: list, mandatory_ids: list) -> dict:
    """Overall status from the registry. PASS iff every mandatory id is present exactly once with state exactly 'PASS'."""
    problems = []
    by_id = {}
    for r in registry:
        by_id.setdefault(r.get("gate_id"), []).append(r)
    states = []
    for gid in mandatory_ids:
        entries = by_id.get(gid, [])
        if not entries:
            problems.append(f"{gid}: missing from registry")
            states.append(FAIL)
            continue
        if len(entries) > 1:
            problems.append(f"{gid}: {len(entries)} registry entries")
            states.append(FAIL)
            continue
        st = entries[0].get("state")
        if st not in STATES:
            problems.append(f"{gid}: unknown state {st!r}")
            states.append(FAIL)
            continue
        states.append(st)
        if st != PASS:
            problems.append(f"{gid}: {st}")
    if not mandatory_ids:
        problems.append("no mandatory gates declared")
        states.append(FAIL)
    if all(s == PASS for s in states) and not problems:
        overall = PASS
    elif FAIL in states:
        overall = FAIL
    elif NOT_EVALUABLE in states:
        overall = NOT_EVALUABLE
    else:
        overall = OVERALL_INCOMPLETE
    return dict(overall_status=overall, n_mandatory=len(mandatory_ids), n_pass=sum(s == PASS for s in states), problems=problems)


def aggregate(specs: list, out_dir: str, extra_gates: Optional[list] = None) -> dict:
    """Evaluate every protocol-derived gate (plus already-evaluated extra gate records) and combine them fail-closed."""
    try:
        cache = {}
        registry = [evaluate_gate(s, out_dir, cache) for s in specs] + list(extra_gates or [])
        mandatory = [s.gate_id for s in specs] + [g["gate_id"] for g in (extra_gates or [])]
        res = combine(registry, mandatory)
        res["registry"] = registry
        return res
    except Exception as e:                                                    # noqa: BLE001
        return dict(overall_status=FAIL, n_mandatory=None, n_pass=None, problems=[f"aggregation exception: {type(e).__name__}: {e}"], registry=[])


def write_registry(path, registry):
    keys = ["gate_id", "state", "reason", "table", "n_expected", "n_rows", "n_missing", "n_duplicate_keys", "n_failed_rows",
            "n_not_evaluable_rows", "n_extra_rows", "first_failures", "description"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in registry:
            w.writerow(r)


# ----------------------------------------------------------------------
# TOST equivalence + power check (plan section 6)
# ----------------------------------------------------------------------
def tost_evaluator(eq_fraction, alpha, under_fraction, n_max):
    from scipy.stats import t as tdist

    def ev(row):
        n, d, se, h = (num(row.get(c)) for c in ("n_k", "d_mean", "se_d", "h"))
        if None in (n, d, se, h) or not all(math.isfinite(x) for x in (n, d, se, h)):
            return "non-numeric control statistics"
        if n != int(n) or n < 2 or n > n_max:
            return f"n_k={n} outside [2, {n_max}]"
        if h <= 0 or se < 0:
            return f"invalid h={h} or se_d={se}"
        eps = eq_fraction * h
        tq = float(tdist.ppf(1.0 - alpha, int(n) - 1))
        if tq * se > under_fraction * eps:
            return f"UNDER-POWERED / NOT DEMONSTRATED: t*SE={tq * se:.3g} > {under_fraction}*eps={under_fraction * eps:.3g}"
        if not (d - tq * se > -eps and d + tq * se < eps):
            return f"TOST not accepted: [{d - tq * se:.3g}, {d + tq * se:.3g}] not inside +-{eps:.3g}"
        return None
    return ev


# ----------------------------------------------------------------------
# B3Y gate specification (from docs/rgf_b3y_protocol.json; read-only regression use)
# ----------------------------------------------------------------------
def b3y_gate_spec(P: dict) -> list:
    L, X = P["limits"], P["xp_gates"]
    S = P["point_sets"]
    pts = {p["id"]: p for p in P["points"]}
    inv = [round(1.0 / d) for d in P["discretization"]["dx_levels"]]                  # [512, 1024, 2048]
    base, mid, fine = inv
    t1, t2, t3v, t3f, buf = S["tier1"], S["tier2_dx2048"], S["tier3_xp_512_validation"], S["tier3_xp_1024_2048"], S["buffer_control"]
    cfg = lambda ids: sorted({pts[i]["config"] for i in ids})
    comp_keys = [(p, base) for p in t1] + [(p, mid) for p in t1] + [(p, fine) for p in t2]
    comp_checks = [Check(c, "le", L["S_abs"]) for c in ("S1_dS", "S2_dS", "S3_dS")] + \
        [Check(c, "le", L["TR_offdiag_abs"]) for c in ("S4_S1_max", "S4_S2_max")] + \
        [Check(f"U_{m}", "le", L["unitarity"]) for m in ("M1", "M2", "M3")] + \
        [Check(f"Q_{m}", "le", L["reciprocity"]) for m in ("M1", "M2", "M3")] + \
        [Check(f"R1_backward_error_{m}", "le", L["sparse_backward_error"]) for m in ("M2", "M3")] + \
        [Check("V_lead_max", "le", L["lead_tail_guard"])] + \
        [Check(c, "true") for c in ("P4_lead_guard_ok", "N1_passed", "P1_rgf_device_bitwise", "P1_M2_samples_bitwise", "P1_M3_samples_bitwise",
                                    "P1_Eperp_bitwise_across_methods")]
    specs = [GateSpec("B3Y.comparison.S1-S6_R1_N1_P1_P4", "Tier 1/2 M1-M3 rows: S1-S6, R1, N1, P1, P4 (no exceedance allowed for PASS)",
                      "per_energy_comparison.csv", ("point", "inv_dx"), expected_keys=comp_keys, checks=comp_checks,
                      status_column="status", status_ok=("OK",), status_not_evaluable=("NOT EVALUABLE",))]
    t3_keys = [(p, base) for p in t3v] + [(p, i) for p in t3f for i in (mid, fine)]
    specs.append(GateSpec("B3Y.tier3.X4_X5_X7_eps", "every M4 solve: X4/X5/X7 gates and eps(Mk, M4) <= limit (X8 at base; adjudication at finer dx)",
                          "tier3_exact_lattice.csv", ("point", "inv_dx"), expected_keys=t3_keys,
                          checks=[Check(c, "true") for c in ("X4_state_passed", "X4_green_passed", "X5_passed", "X7_passed", "gates_passed")] +
                          [Check("X7_error_indicator", "le", X["error_indicator_max"]), Check("X5_fl_vs_state_max_abs", "le", X["fisher_lee_vs_wavefunction_forms_abs"])] +
                          [Check(f"eps_{m}", "le", L["S_abs"]) for m in ("M1", "M2", "M3")]))
    dec_keys = [(p, i, m) for p in t3f for i in (mid, fine) for m in ("M2", "M3")]
    specs.append(GateSpec("B3Y.decomposition.evaluable", "every protocol-required M2/M3 decomposition row is evaluable (gates, mirrors, finite terms)",
                          "error_decomposition.csv", ("point", "inv_dx", "method"), expected_keys=dec_keys,
                          status_column="status", status_ok=("", "OK"), status_not_evaluable=("NOT EVALUABLE",),
                          checks=[Check(c, "finite") for c in ("total", "bulk_rounding", "lead_data", "solver", "identity_residual_three_term", "e_ind_combined")] +
                          [Check(c, "true") for c in ("gate_fl_passed", "gate_hybrid_passed", "mirror_source_bitwise", "mirror_extraction_bitwise",
                                                      "M2_recomputed_bitwise", "M3_recomputed_bitwise")]))
    xg = [("G0", [(p, "") for p in t1], ("case", "inv_dx")), ("X1", [(p, i) for p in dict.fromkeys(t3v + t3f) for i in inv], ("case", "inv_dx")),
          ("X2", [(p, i) for p in t3f for i in inv], ("case", "inv_dx")), ("X3", [(p, i) for p in t3f for i in inv], ("case", "inv_dx")),
          ("X6", [(p, i) for p in t3f for i in (base, mid)], ("case", "inv_dx"))]
    for g, keys, kc in xg:
        specs.append(GateSpec(f"B3Y.xp_gate.{g}", f"gate {g}: every protocol-required row passed", "xp_gates.csv", kc, expected_keys=keys,
                              row_filter={"gate": g}, checks=[Check("passed", "true")]))
    for g in ("X0", "X8"):
        specs.append(GateSpec(f"B3Y.xp_gate.{g}", f"gate {g}: exactly one row, passed", "xp_gates.csv", ("gate",), expected_count=1,
                              row_filter={"gate": g}, checks=[Check("passed", "true")]))
    specs.append(GateSpec("B3Y.K1_buffer_control", "K1 buffer control rows", "buffer_control.csv", ("point",), expected_keys=[(p,) for p in buf],
                          checks=[Check("K1_passed_dT1", "true"), Check("S_within_geometry_passed", "true"),
                                  Check("abs_dT1_M1", "le", L["dT1_buffer"]), Check("abs_dT1_M3", "le", L["dT1_buffer"])]))
    specs.append(GateSpec("B3Y.H1_solver_spread", "H1 T1 spread rows", "solver_spread.csv", ("point", "inv_dx"), expected_keys=[(p, i) for p in t2 for i in inv],
                          checks=[Check("passed", "true"), Check("T1_spread", "le", L["T1_method_spread_for_solver_independence"])]))
    specs.append(GateSpec("B3Y.O1_rows_present", "O1 convergence rows present (classification only; completeness gate)", "refinement_convergence.csv",
                          ("point", "method"), expected_keys=[(p, m) for p in t2 for m in ("M1", "M2", "M3", "M4")], checks=[Check("observed_order", "finite")]))
    dev = [(c, base, 0.2) for c in cfg(t1)] + [(c, mid, 0.2) for c in cfg(t1)] + [(c, fine, 0.2) for c in cfg(t2)] + [(c, mid, 0.4) for c in cfg(buf)]
    specs.append(GateSpec("B3Y.P1_P4_devices", "potential identity and lead guard of every predeclared device", "potential_identity.csv",
                          ("config", "inv_dx", "buffer_requested"), expected_keys=dev, allow_duplicates=True,
                          row_filter={},
                          checks=[]))
    # P1/P4 checks only on device rows (P2 rows carry an inv_dx like '1024_vs_512'); a dedicated evaluator keeps P2 rows out of P1 checks
    def p1_eval(row):
        if "_vs_" in str(row.get("inv_dx", "")):
            return None
        for c in ("P1_rgf_device_vs_direct_bitwise", "P4_lead_guard_ok"):
            if not is_true(row.get(c)):
                return f"{c}={row.get(c)!r}"
        v = num(row.get("V_lead_max"))
        return None if (v is not None and v <= L["lead_tail_guard"]) else f"V_lead_max={row.get('V_lead_max')!r}"
    specs[-1].row_evaluator = p1_eval
    specs.append(GateSpec("B3Y.P2_cross_resolution", "coincident-slice potential identity across resolutions", "potential_identity.csv",
                          ("config", "inv_dx"), expected_keys=[(c, f"{i}_vs_{base}") for c in cfg(t2) for i in (mid, fine)],
                          row_filter={"buffer_requested": ""}, checks=[Check("P2_potential_bitwise", "true"), Check("coincident_x_bitwise", "true")]))
    return specs


def aggregate_b3y(out_dir: str, protocol_path: Optional[str] = None) -> dict:
    """Read-only re-aggregation of existing B3Y outputs (regression of the fail-closed logic). Writes nothing."""
    protocol_path = protocol_path or os.path.join(HERE, "docs", "rgf_b3y_protocol.json")
    with open(protocol_path) as f:
        P = json.load(f)
    nested = nested_status_gate("B3Y.summary.nested_statuses", os.path.join(out_dir, "summary.json"), "criteria")
    return aggregate(b3y_gate_spec(P), out_dir, extra_gates=[nested])


# ----------------------------------------------------------------------
# physical-campaign gate specification (from the DRAFT gate file)
# ----------------------------------------------------------------------
def _cells(D, which):
    if which == "all":
        return [(w, z) for w in D["grid"]["W"] for z in D["grid"]["zeta"]]
    if which == "corners":
        return [tuple(c) for c in D["grid"]["corner_cells"]]
    return [tuple(c) for c in which]


def cell_id(w, z):
    return f"W{float(w)!r}_z{float(z)!r}"


def physical_gate_spec(D: dict, n_seeds: int) -> list:
    S, L = D["seeds"], D["limits"]
    if n_seeds is None or not (S["N_min"] <= int(n_seeds) <= S["N_max"]) or int(n_seeds) != n_seeds:
        raise ValueError(f"N = {n_seeds!r} is not a protocol-admissible confirmatory ensemble size [{S['N_min']}, {S['N_max']}]")
    N = int(n_seeds)
    seeds = list(range(S["confirmatory_first"], S["confirmatory_first"] + N))
    E = list(range(D["energies"]["n_energies"]))
    cells = _cells(D, "all")
    C = D["controls"]
    rs = D["row_status_vocabulary"]
    specs = [GateSpec("PHYS.stages.completed", "every predeclared stage completed (no early termination)", "stages.csv", ("stage",),
                      expected_keys=[(s,) for s in D["stages"]], status_column="status", status_ok=("COMPLETED",), extra_rows_allowed=False)]
    specs.append(GateSpec("PHYS.C0.realizations", "C0 per-realization gates: unitarity, reciprocity, channel count, lead tail, V hash",
                          "realizations.csv", ("W", "zeta", "seed", "energy_index"),
                          expected_keys=[(w, z, s, e) for (w, z) in cells for s in seeds for e in E],
                          status_column="status", status_ok=tuple(rs["ok"]), status_not_evaluable=tuple(rs["not_evaluable"]),
                          checks=[Check("unitarity", "le", L["unitarity"]), Check("reciprocity", "le", L["reciprocity"]),
                                  Check("lead_tail_max", "le", L["lead_tail_guard"]), Check("n_open", "eqcol", "n_open_analytic"),
                                  Check("v_hash", "nonempty"), Check("tau", "finite"), Check("T1", "finite")]))
    specs.append(GateSpec("PHYS.C0.reproducibility", "C0 bitwise-reproducible rerun of a predeclared-size random subset", "reproducibility.csv",
                          ("W", "zeta", "seed"), expected_count=math.ceil(L["reproducibility_subset_fraction"] * len(cells) * N),
                          checks=[Check("bitwise_equal", "true")]))
    c2 = C["C2"]
    specs.append(GateSpec("PHYS.C2.solver_cross_check", "C2 M3 vs M1 max|dS| (unadjudicated or adjudicated exceedance is never PASS)",
                          "solver_cross_check.csv", ("W", "zeta", "seed", "energy_index"),
                          expected_keys=[(w, z, s, e) for (w, z) in _cells(D, c2["cells"]) for s in seeds[:c2["n_seeds"]] for e in E],
                          checks=[Check("max_abs_dS", "le", L["solver_cross_check_S_abs"])]))
    specs.append(GateSpec("PHYS.C2.adjudication_evaluable", "every adjudication row (if any) is evaluable", "adjudication.csv",
                          ("W", "zeta", "seed", "energy_index"), expected_keys=[], status_column="status", status_ok=("ADJUDICATED",),
                          status_not_evaluable=("NOT EVALUABLE",)))
    c5 = C["C5"]
    specs.append(GateSpec("PHYS.C5.b3y_reproduction", "C5 seed 0 reproduces B3Y", "b3y_reproduction.csv", ("zeta", "energy_index"),
                          expected_keys=[(z, e) for z in c5["zeta"] for e in E], checks=[Check("abs_dT1", "le", L["b3y_reproduction_abs_dT1"])]))
    specs.append(GateSpec("PHYS.C6.clean_guide", "C6 clean guide T_nn = 1", "clean_guide.csv", ("energy_index",), expected_keys=[(e,) for e in E],
                          checks=[Check("max_abs_Tnn_minus_1", "le", L["clean_guide_Tnn_abs"])]))
    tost = tost_evaluator(L["equivalence_fraction_of_h"], L["tost_alpha"], L["underpowered_fraction_of_eps"], L["control_seed_extension_max"])
    contrasts = D["contrasts"]
    for cid in ("C1", "C1b", "C3", "C4"):
        c = C[cid]
        cs = set(_cells(D, c["cells"]))
        keys = []
        for o in c["observables"]:
            keys += [(cid, o, "cell", cell_id(w, z)) for (w, z) in sorted(cs)]
            keys += [(cid, o, "contrast", k["id"]) for k in contrasts if all(tuple(x) in cs for x in k["cells"])]
        specs.append(GateSpec(f"PHYS.{cid}.equivalence", f"{cid} {c['what']}: TOST equivalence at alpha = {L['tost_alpha']} with power check, per cell and covered contrast",
                              "controls.csv", ("control_id", "observable", "unit_type", "unit_id"), expected_keys=keys,
                              row_filter={"control_id": cid}, row_evaluator=tost))
    c7 = C["C7"]
    keys = [(c7["observable"], "cell", cell_id(w, z)) for (w, z) in cells] + [(c7["observable"], "contrast", k["id"]) for k in contrasts]
    specs.append(GateSpec("PHYS.C7.precision", "C7 achieved 95 % CI half-width <= delta (primary observable only)", "precision.csv",
                          ("observable", "unit_type", "unit_id"), expected_keys=keys, row_filter={"observable": c7["observable"]},
                          checks=[Check("half_width", "le", L["delta"])]))
    dec = tuple(D["statistics"]["decisions"])

    def test_eval(row):
        p = num(row.get("p_holm_adj_from_phi"))
        if p is None or not (0.0 <= p <= 1.0):
            return f"p_holm_adj_from_phi={row.get('p_holm_adj_from_phi')!r} not in [0, 1]"
        d = row.get("decision")
        if d == dec[0] and not p <= L["holm_alpha"]:
            return "REJECTED with adjusted p above alpha"
        if d == dec[1] and p <= L["holm_alpha"]:
            return "NOT REJECTED with adjusted p <= alpha"
        return None
    specs.append(GateSpec("PHYS.statistics.tests_complete", "all 68 confirmatory tests present with a valid decision", "tests.csv",
                          ("observable", "contrast_id"), expected_keys=[(o, k["id"]) for o in D["statistics"]["observables"] for k in contrasts],
                          status_column="decision", status_ok=dec, row_evaluator=test_eval))
    if len(specs[-1].expected_keys) != D["statistics"]["n_tests"]:
        raise ValueError("draft: number of tests inconsistent with the contrast list")
    return specs


def load_draft(path=None):
    with open(path or os.path.join(HERE, "docs", "physical_campaign_gates_draft.json")) as f:
        return json.load(f)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="read-only fail-closed re-aggregation")
    ap.add_argument("--b3y", required=True, help="existing B3Y result directory (read only)")
    a = ap.parse_args()
    r = aggregate_b3y(a.b3y)
    for g in r["registry"]:
        print(f"{g['state']:14s} {g['gate_id']:40s} expected={g['n_expected']} rows={g['n_rows']} {g['reason']}")
    print("OVERALL:", r["overall_status"], r["problems"])
