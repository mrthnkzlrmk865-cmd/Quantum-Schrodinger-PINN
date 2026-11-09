r"""
Fail-closed gate registry and overall-status aggregation for validation / physical campaigns.

Design (docs/PHYSICAL_CAMPAIGN_PLAN.md section 8):
  * The registry of MANDATORY gates is built from the predeclared protocol (exact row-key sets, thresholds), never from whatever rows a
    run happened to produce.
  * Every gate has exactly one state in STATES = (PASS, FAIL, NOT RUN, NOT EVALUABLE). Overall PASS iff every mandatory gate is
    explicitly PASS. A missing gate, an unknown state, an exception, a missing table, a missing / duplicated / UNEXPECTED row, a
    malformed CSV row, an unknown row status, or a stage that did not complete all prevent PASS.
  * Row-key sets are exact: the observed keys must equal the protocol-derived keys (optional rows only where the protocol defines them
    explicitly via `allowed_keys`; repeated rows only where the protocol defines their multiplicity). Rows that no gate claims fail the
    per-table coverage gate. Nothing is ever filtered out of the output tables.
  * Statuses inside nested upstream summaries are flattened recursively.

Gate-specification builders:
  b3y_gate_spec(protocol)                 from docs/rgf_b3y_protocol.json (read-only use: regression against existing B3Y results)
  physical_gate_spec(draft, n_seeds)      from docs/physical_campaign_gates_draft.json (DRAFT; not frozen)

This module does not import, and does not depend on, validate_rgf_b3y.py.
Command line (read-only):  python campaign_status.py --b3y results/rgf_b3y_validation
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from collections import Counter
from fractions import Fraction
from dataclasses import dataclass, field
from typing import Callable, Optional

PASS, FAIL, NOT_RUN, NOT_EVALUABLE = "PASS", "FAIL", "NOT RUN", "NOT EVALUABLE"
STATES = (PASS, FAIL, NOT_RUN, NOT_EVALUABLE)
OVERALL_INCOMPLETE = "INCOMPLETE"
OVERALL_STATES = (PASS, FAIL, NOT_EVALUABLE, OVERALL_INCOMPLETE)
HERE = os.path.dirname(os.path.abspath(__file__))

# Decision registries. Each maps an APPROVED rule name to this module's own validator; the independent verifier keeps separate ones.
#   SELECTION_ALGORITHMS: C0 subset selection, name -> function(ensemble_units, n, seed, spec) -> list of selected units
#   RNG_RULES:            bootstrap stream rule, name -> function(meta, spec, N) -> list of failure reasons
# They are filled at the end of this module with the rules approved in R7 (docs/CAMPAIGN_STATUS_CONTRACT.md section 1 and 5).
SELECTION_ALGORITHMS: dict = {}
RNG_RULES: dict = {}
# C2 adjudication criteria. EMPTY on purpose: no numerical resolution criterion is approved (draft: C2_adjudication_criteria = null).
# Name -> function(adjudication_row, criteria_spec) -> None | reason.
C2_ADJUDICATION_CRITERIA: dict = {}

# Numerical verification tolerance for REPORTED floating-point statistics (section 5.3): a and b agree iff both are finite and
#   |a - b| <= VERIFY_ATOL + VERIFY_RTOL * max(|a|, |b|).
# Justification (docs/CAMPAIGN_STATUS_CONTRACT.md, section 2): the largest legitimate implementation difference measured is 2.5e-11
# relative (mc_se at K = B-1, B = 1e6, cancellation in 1 - p_hat); scipy Clopper-Pearson endpoints agree with exact binomial tails to
# <= 7.6e-14 relative. ATOL = 0 because every compared quantity is either structurally exact (0 at K = 0, 1 at K = B; checked exactly)
# or >= 3.7e-11 (smallest nonzero CP lower endpoint), where any positive ATOL would only weaken the test. Decisions never depend on
# this tolerance: Holm sets are recomputed from the integer counts.
VERIFY_RTOL = 1e-9
VERIFY_ATOL = 0.0


def approx_equal(a, b, rtol=VERIFY_RTOL, atol=VERIFY_ATOL) -> bool:
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return False
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= atol + rtol * max(abs(a), abs(b))


# ----------------------------------------------------------------------
# value normalisation
# ----------------------------------------------------------------------
def norm(v):
    """Canonical comparison value: numbers compare as floats ('25', '25.0', 25 are equal), everything else as stripped text.
    None (a missing CSV field) is kept distinct from every real value."""
    if v is None:
        return None
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
    """Strict float parse; None for None, booleans, empty or non-numeric text."""
    if isinstance(v, bool) or v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def is_true(v):
    return v is True or (isinstance(v, str) and v.strip() == "True")


def is_nonempty(v):
    """None, '' and whitespace-only text are empty (str(None) == 'None' must NOT count as content)."""
    return v is not None and isinstance(v, (str, int, float)) and not isinstance(v, bool) and str(v).strip() != ""


# ----------------------------------------------------------------------
# checks and gate specifications
# ----------------------------------------------------------------------
@dataclass
class Check:
    column: str
    op: str                  # 'le' (value <= limit), 'eq' (numeric equality), 'true', 'finite', 'nonempty', 'eqcol' (== another column)
    limit: object = None

    def evaluate(self, row) -> Optional[str]:
        """None if the check passes, otherwise a reason. Missing columns, None and unparsable values fail (fail-closed)."""
        if self.column not in row:
            return f"column {self.column!r} missing"
        v = row[self.column]
        if self.op == "true":
            return None if is_true(v) else f"{self.column}={v!r} is not True"
        if self.op == "nonempty":
            return None if is_nonempty(v) else f"{self.column}={v!r} empty"
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
    expected_keys: Optional[list] = None               # protocol-derived REQUIRED keys, each exactly once
    allowed_keys: Optional[list] = None                # protocol-defined OPTIONAL keys (at most once each); anything else is unexpected
    expected_multiplicity: Optional[dict] = None       # protocol-derived key -> exact number of rows (tables that repeat a record)
    expected_keys_fn: Optional[Callable] = None        # (load) -> required keys, for protocol rules that depend on another table
    allowed_keys_fn: Optional[Callable] = None         # (load) -> optional keys whose permission depends on another table
    row_filter: dict = field(default_factory=dict)     # column -> value; selects the rows owned by this gate
    row_select: Optional[Callable] = None              # additional ownership predicate
    checks: list = field(default_factory=list)
    status_column: Optional[str] = None
    status_ok: tuple = ()
    status_not_evaluable: tuple = ()
    row_evaluator: Optional[Callable] = None           # row -> None | reason
    table_evaluator: Optional[Callable] = None         # rows -> list of reasons (cross-row rules, e.g. Holm sets)
    not_evaluable_reason: Optional[str] = None         # a missing protocol decision: the gate cannot be PASS until it is made
    not_evaluable_fn: Optional[Callable] = None        # (load) -> reason | None (a missing decision that matters only for some data)

    def owns(self, row) -> bool:
        if any(norm(row.get(c)) != norm(v) for c, v in self.row_filter.items()):
            return False
        return self.row_select(row) if self.row_select is not None else True


def _row_key(row, cols):
    return tuple(norm(row.get(c)) for c in cols)


def read_table(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def malformed(row) -> Optional[str]:
    """csv.DictReader stores surplus fields under the key None and missing trailing fields as the value None."""
    if None in row:
        return "row has more fields than the header"
    if any(v is None for v in row.values()):
        return "row has fewer fields than the header"
    return None


def _empty_record(spec_id, description, table):
    return dict(gate_id=spec_id, description=description, table=table, state=FAIL, reason="", n_expected=None, n_rows=0, n_missing=0,
                n_duplicate_keys=0, n_failed_rows=0, n_not_evaluable_rows=0, n_extra_rows=0, first_failures="")


def _load(out_dir, table, cache):
    if table not in cache:
        p = os.path.join(out_dir, table)
        cache[table] = read_table(p) if os.path.exists(p) else None
    return cache[table]


def evaluate_gate(spec: GateSpec, out_dir: str, cache: dict) -> dict:
    """Evaluate one mandatory gate. Never raises: any exception is a FAIL of that gate."""
    rec = _empty_record(spec.gate_id, spec.description, spec.table)
    try:
        table = _load(out_dir, spec.table, cache)
        if table is None:
            rec.update(state=NOT_RUN, reason=f"table {spec.table} missing")
            return rec
        rows = [r for r in table if spec.owns(r)]
        rec["n_rows"] = len(rows)
        reasons = []
        counts = Counter(_row_key(r, spec.key_columns) for r in rows)
        if spec.expected_multiplicity is not None:
            exp = Counter({tuple(norm(x) for x in k): n for k, n in spec.expected_multiplicity.items()})
            rec["n_expected"] = sum(exp.values())
            missing = [k for k in exp if counts.get(k, 0) < exp[k]]
            dup = [k for k in exp if counts.get(k, 0) > exp[k]]
            extra = [k for k in counts if k not in exp]
        else:
            load = lambda t: _load(out_dir, t, cache)
            req = spec.expected_keys_fn(load) if spec.expected_keys_fn is not None else spec.expected_keys
            if req is None:
                reasons.append("gate has no protocol-derived key set")
                req = []
            req = [tuple(norm(x) for x in k) for k in req]
            opt = [tuple(norm(x) for x in k) for k in (spec.allowed_keys or [])]
            if spec.allowed_keys_fn is not None:
                opt += [tuple(norm(x) for x in k) for k in spec.allowed_keys_fn(load)]
            if len(set(req)) != len(req) or set(req) & set(opt):
                reasons.append("protocol-derived key set is not unique")
            rec["n_expected"] = len(req)
            allowed = set(req) | set(opt)
            missing = [k for k in req if k not in counts]
            dup = [k for k, n in counts.items() if n > 1]
            extra = [k for k in counts if k not in allowed]
        rec.update(n_missing=len(missing), n_duplicate_keys=len(dup), n_extra_rows=sum(counts[k] for k in extra))
        if missing:
            reasons.append(f"{len(missing)} expected rows missing (e.g. {missing[0]})")
        if dup:
            reasons.append(f"{len(dup)} keys with too many rows (e.g. {dup[0]})")
        if extra:
            reasons.append(f"{len(extra)} unexpected keys (e.g. {extra[0]})")
        failed, ne = [], 0
        for r in rows:
            bad = malformed(r)
            if bad:
                failed.append(f"{_row_key(r, spec.key_columns)}: {bad}")
                continue
            if spec.status_column is not None:
                st = r[spec.status_column] if spec.status_column in r else ""      # column absent from the table: only '' (if declared ok) passes
                st = st.strip() if isinstance(st, str) else st
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
        if spec.table_evaluator is not None and not failed and not reasons:
            failed += list(spec.table_evaluator(rows))
        rec["n_failed_rows"], rec["n_not_evaluable_rows"] = len(failed), ne
        if failed:
            reasons.append(f"{len(failed)} row/table failures")
            rec["first_failures"] = " | ".join(failed[:3])
        if reasons:
            rec.update(state=FAIL, reason="; ".join(reasons))
        elif ne:
            rec.update(state=NOT_EVALUABLE, reason=f"{ne} NOT EVALUABLE rows")
        elif spec.not_evaluable_reason:
            rec.update(state=NOT_EVALUABLE, reason=spec.not_evaluable_reason)
        elif spec.not_evaluable_fn is not None and spec.not_evaluable_fn(lambda t: _load(out_dir, t, cache)):
            rec.update(state=NOT_EVALUABLE, reason=spec.not_evaluable_fn(lambda t: _load(out_dir, t, cache)))
        else:
            rec.update(state=PASS, reason="observed key set equals the protocol key set; all rows pass")
    except Exception as e:                                                    # noqa: BLE001  (fail closed)
        rec.update(state=FAIL, reason=f"exception during evaluation: {type(e).__name__}: {e}")
    return rec


def coverage_gates(specs: list, out_dir: str, cache: dict) -> list:
    """One gate per table: every row of the table must be owned by at least one gate (no unclaimed rows)."""
    out = []
    for t in dict.fromkeys(s.table for s in specs):
        rec = _empty_record(f"TABLE.{t}.coverage", f"every row of {t} is claimed by a protocol gate", t)
        try:
            table = _load(out_dir, t, cache)
            if table is None:
                rec.update(state=NOT_RUN, reason="table missing")
            else:
                owners = [s for s in specs if s.table == t]
                orphan = [r for r in table if not any(s.owns(r) for s in owners)]
                rec.update(n_rows=len(table), n_extra_rows=len(orphan))
                if orphan:
                    rec.update(state=FAIL, reason=f"{len(orphan)} rows not claimed by any gate", first_failures=str(dict(orphan[0]))[:200])
                else:
                    rec.update(state=PASS, reason="all rows claimed")
        except Exception as e:                                                # noqa: BLE001
            rec.update(state=FAIL, reason=f"exception: {type(e).__name__}: {e}")
        out.append(rec)
    return out


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
    """Gate on an upstream summary: every nested status PASS (or an explicitly informational label), every nested 'passed' True,
    and the early-termination fields null."""
    rec = _empty_record(gate_id, f"nested statuses of {os.path.basename(summary_path)}:{subtree}", os.path.basename(summary_path))
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
        flat = flatten_statuses(s[subtree], subtree) if subtree in s else []
        if subtree not in s:
            bad.append(f"{subtree} missing")
        elif not flat:
            bad.append("no statuses found")
        rec["n_rows"] = len(flat)
        for p, v in flat:
            if p.endswith("/passed"):
                if v is not True:
                    bad.append(f"{p} = {v!r}")
            elif not (v == PASS or v in informational):
                bad.append(f"{p} = {v!r}")
        rec["n_failed_rows"] = len(bad)
        if bad:
            ne_only = all(b.endswith("= 'NOT EVALUABLE'") for b in bad)
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
        if len(entries) != 1:
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
    """Evaluate every protocol-derived gate, the per-table coverage gates and any extra gate records; combine fail-closed."""
    try:
        cache = {}
        registry = [evaluate_gate(s, out_dir, cache) for s in specs] + coverage_gates(specs, out_dir, cache) + list(extra_gates or [])
        mandatory = [g["gate_id"] for g in registry]
        if len(set(mandatory)) != len(mandatory):
            raise ValueError("duplicate gate ids in the specification")
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
# section 5.3 statistical output contract (validation only; the bootstrap itself is not implemented here)
# ----------------------------------------------------------------------
def clopper_pearson(K, B, alpha):
    from scipy.stats import beta
    lo = 0.0 if K == 0 else float(beta.ppf(alpha / 2.0, K, B - K + 1))
    hi = 1.0 if K == B else float(beta.ppf(1.0 - alpha / 2.0, K + 1, B - K))
    return lo, hi


def holm_set(p: dict, alpha: float) -> set:
    """Holm step-down rejection set for {test: p}."""
    order = sorted(p, key=lambda k: p[k])
    m, out = len(order), set()
    for i, k in enumerate(order):
        if p[k] <= alpha / (m - i):
            out.add(k)
        else:
            break
    return out


def statistics_table_evaluator(D: dict, test_order: list):
    st = D["statistics"]
    alpha, B1, B2 = st["holm_alpha"], st["B_stage1"], st["B_stage2"]
    a1, a2 = st["mc_alpha_per_test_stage1"], st["mc_alpha_per_test_stage2"]
    REJ, NREJ, AMB = st["decisions"]
    rel = approx_equal

    def ev(rows):
        bad = []
        V = {}
        for r in rows:
            key = (r["observable"], r["contrast_id"])
            try:
                K1, K, B, stage = (int(float(r[c])) for c in ("K_stage1", "K", "B", "stage"))
                assert float(r["K"]) == K and float(r["K_stage1"]) == K1
                f = {c: float(r[c]) for c in ("theta_hat", "ci95_sampling_lo", "ci95_sampling_hi", "p_hat", "mc_se", "mc_ci_sim99_lo",
                                               "mc_ci_sim99_hi", "p_lo_final", "p_hi_final", "holm_rank", "holm_threshold", "p_holm_adj_from_phi")}
            except (KeyError, TypeError, ValueError, AssertionError) as e:
                bad.append(f"{key}: unparsable statistics ({type(e).__name__})")
                continue
            if not all(math.isfinite(v) for v in f.values()):
                bad.append(f"{key}: non-finite value")
                continue
            if not f["ci95_sampling_lo"] <= f["ci95_sampling_hi"]:
                bad.append(f"{key}: sampling CI reversed")
            if (stage, B) not in ((1, B1), (2, B2)) or not (0 <= K1 <= B1) or not (0 <= K <= B) or (stage == 1 and K != K1):
                bad.append(f"{key}: inconsistent stage/B/K")
                continue
            ph = (K + 1) / (B + 1)
            lo, hi = clopper_pearson(K, B, a1 if stage == 1 else a2)
            if not (rel(f["p_hat"], ph) and rel(f["mc_se"], math.sqrt(ph * (1 - ph) / B)) and rel(f["mc_ci_sim99_lo"], lo)
                    and rel(f["mc_ci_sim99_hi"], hi) and rel(f["p_lo_final"], lo) and rel(f["p_hi_final"], max(hi, ph))):
                bad.append(f"{key}: p_hat / MC SE / simultaneous MC interval / endpoint vectors inconsistent with K, B")
                continue
            if (K == 0 and (f["mc_ci_sim99_lo"] != 0.0 or f["p_lo_final"] != 0.0)) or (K == B and (f["mc_ci_sim99_hi"] != 1.0 or f["p_hi_final"] != 1.0)):
                bad.append(f"{key}: structurally exact endpoint (0 at K = 0, 1 at K = B) not exact")
                continue
            if not (0.0 <= f["p_lo_final"] <= f["p_hat"] <= f["p_hi_final"] <= 1.0):
                bad.append(f"{key}: endpoint ordering 0 <= p_lo <= p_hat <= p_hi <= 1 violated")
                continue
            lo1, hi1 = clopper_pearson(K1, B1, a1)
            # Holm sets, ranks and adjusted p are computed from the values implied by the INTEGER counts, never from the reported floats
            V[key] = dict(stage=stage, plo1=lo1, phi1=max(hi1, (K1 + 1) / (B1 + 1)), plo=lo, phi=max(hi, ph), f=f, d=r["decision"])
        if bad:
            return bad
        if set(V) != set(test_order):
            return ["statistics rows do not match the 68-test family"]
        A1 = holm_set({k: V[k]["plo1"] for k in V}, alpha) - holm_set({k: V[k]["phi1"] for k in V}, alpha)
        for k in V:
            if (V[k]["stage"] == 2) != (k in A1):
                bad.append(f"{k}: stage {V[k]['stage']} but Stage-1 ambiguity membership is {k in A1}")
        R = holm_set({k: V[k]["phi"] for k in V}, alpha)
        A = holm_set({k: V[k]["plo"] for k in V}, alpha) - R
        idx = {k: i for i, k in enumerate(test_order)}
        order = sorted(V, key=lambda k: (V[k]["phi"], idx[k]))
        m, run = len(order), 0.0
        for i, k in enumerate(order):
            run = max(run, min(1.0, (m - i) * V[k]["phi"]))
            f = V[k]["f"]
            if f["holm_rank"] != i + 1 or not rel(f["holm_threshold"], alpha / (m - i)) or not rel(f["p_holm_adj_from_phi"], run):
                bad.append(f"{k}: Holm rank/threshold/adjusted p inconsistent with p_hi_final")
            want = REJ if k in R else AMB if k in A else NREJ
            if V[k]["d"] != want:
                bad.append(f"{k}: decision {V[k]['d']!r}, Holm sets require {want!r}")
            if (k in R) != (run <= alpha):
                bad.append(f"{k}: adjusted p and R_claim disagree")
        return bad
    return ev


# ----------------------------------------------------------------------
# bootstrap RNG stream metadata (section 5.3 (a), (f) step 3; docs/CAMPAIGN_STATUS_CONTRACT.md section 1)
# ----------------------------------------------------------------------
RNG_COLUMNS = ("bootstrap_rng_rule", "bootstrap_rng_stream_stage1", "bootstrap_seed_stage1", "bootstrap_rng_stream_stage2", "bootstrap_seed_stage2",
               "bootstrap_bit_generator", "bootstrap_numpy_version", "bootstrap_resampling_api", "bootstrap_probe_stage1", "bootstrap_probe_stage2")
STAGE_INVARIANT = ("bootstrap_rng_rule", "bootstrap_bit_generator", "bootstrap_numpy_version", "bootstrap_resampling_api")


def _seed_ok(v) -> bool:
    """A recorded seed / entropy must be a non-negative integer written in decimal digits (no sign, no decimal point, not blank)."""
    return isinstance(v, str) and v.strip().isdigit()


def numpy_version() -> str:
    import numpy
    return numpy.__version__


def bootstrap_rng_evaluator(D: dict, N: int):
    """Checks that hold for ANY stream rule (one shared Stage-1 stream for the paired replicates of all 68 tests; one shared Stage-2 stream for
    all Stage-2 tests; distinct stream identifiers; Stage-2 fields empty exactly for Stage-1 rows; stage-invariant metadata identical in
    every row), then the approved rule's own validator."""
    spec = D.get("bootstrap_rng")

    def ev(rows):
        bad = []
        inv, s1, s2 = set(), set(), set()
        for r in rows:
            key = (r.get("observable"), r.get("contrast_id"))
            if any(c not in r for c in RNG_COLUMNS):
                bad.append(f"{key}: RNG metadata columns missing")
                continue
            stage = norm(r.get("stage"))
            if not (all(is_nonempty(r[c]) for c in STAGE_INVARIANT) and is_nonempty(r["bootstrap_rng_stream_stage1"])
                    and _seed_ok(r["bootstrap_seed_stage1"]) and is_nonempty(r["bootstrap_probe_stage1"])):
                bad.append(f"{key}: Stage-1 RNG metadata missing or malformed")
                continue
            inv.add(tuple(r[c].strip() for c in STAGE_INVARIANT))
            s1.add((r["bootstrap_rng_stream_stage1"].strip(), int(r["bootstrap_seed_stage1"]), r["bootstrap_probe_stage1"].strip()))
            has2 = any(is_nonempty(r[c]) for c in ("bootstrap_rng_stream_stage2", "bootstrap_seed_stage2", "bootstrap_probe_stage2"))
            if stage == 1.0:
                if has2:
                    bad.append(f"{key}: Stage-1 test carries Stage-2 RNG metadata")
            elif stage == 2.0:
                if not (is_nonempty(r["bootstrap_rng_stream_stage2"]) and _seed_ok(r["bootstrap_seed_stage2"]) and is_nonempty(r["bootstrap_probe_stage2"])):
                    bad.append(f"{key}: Stage-2 RNG metadata missing or malformed")
                    continue
                s2.add((r["bootstrap_rng_stream_stage2"].strip(), int(r["bootstrap_seed_stage2"]), r["bootstrap_probe_stage2"].strip()))
            else:
                bad.append(f"{key}: stage {r.get('stage')!r}")
        if bad:
            return bad
        if len(inv) != 1:
            bad.append(f"{len(inv)} different (rule, bit generator, numpy version, resampling API) records")
        if len(s1) != 1:
            bad.append(f"{len(s1)} different Stage-1 streams recorded (the paired Stage-1 replicates are shared by all tests)")
        if len(s2) > 1:
            bad.append(f"{len(s2)} different Stage-2 streams recorded (the Stage-2 replicates are shared by all Stage-2 tests)")
        if s1 and s2 and next(iter(s1))[0] == next(iter(s2))[0]:
            bad.append("Stage-2 stream identifier equals the Stage-1 stream identifier")
        if bad or spec is None:
            return bad
        rule, bitgen, npv, api = next(iter(inv))
        if rule != str(spec.get("rule")):
            return [f"recorded stream rule {rule!r} is not the approved rule {spec.get('rule')!r}"]
        fn = RNG_RULES.get(spec.get("rule"))
        if fn is None or spec.get("root_entropy") is None:
            return bad                                                   # structure only; the gate is NOT EVALUABLE (decision missing)
        meta = dict(rule=rule, bit_generator=bitgen, numpy_version=npv, resampling_api=api, stage1=next(iter(s1)), stage2=next(iter(s2)) if s2 else None)
        return list(fn(meta, spec, N))
    return ev


def bootstrap_rng_gate(D: dict, order: list, N: int) -> GateSpec:
    spec = D.get("bootstrap_rng")
    if spec is None:
        reason = "PROTOCOL DECISION REQUIRED: bootstrap RNG rule not specified (independence by design not established)"
    elif spec.get("rule") not in RNG_RULES:
        reason = f"no validator implemented for the bootstrap RNG rule {spec.get('rule')!r}"
    elif spec.get("root_entropy") is None:
        reason = "PROTOCOL DECISION REQUIRED: the root SeedSequence entropy of the approved rule has not been recorded"
    elif not _seed_ok(str(spec.get("root_entropy"))):
        raise ValueError(f"protocol root_entropy {spec.get('root_entropy')!r} is not a non-negative decimal integer")
    else:
        reason = None

    def numpy_pin(load):
        if spec is not None and spec.get("numpy_version") and numpy_version() != spec["numpy_version"]:
            return f"installed numpy {numpy_version()} differs from the pinned {spec['numpy_version']}: streams cannot be reproduced bit for bit here"
        return None
    return GateSpec("PHYS.statistics.bootstrap_rng", "bootstrap RNG: approved SeedSequence/PCG64 spawn rule, Stage-1 and Stage-2 streams recorded and reproducible",
                    "tests.csv", ("observable", "contrast_id"), expected_keys=order, table_evaluator=bootstrap_rng_evaluator(D, N),
                    not_evaluable_reason=reason, not_evaluable_fn=numpy_pin)


# ----- approved rule (R7, option A): SeedSequence(entropy=E) root, PCG64, Stage 1 = child spawn_key (0,), Stage 2 = child spawn_key (1,) -----
RULE_A = "SEEDSEQ-PCG64-SPAWN-A"


def rule_a_stream_id(E: int, stage: int, numpy_ver: str) -> str:
    return f"{RULE_A}|entropy={int(E)}|spawn_key=({stage - 1},)|bitgen=PCG64|numpy={numpy_ver}"


def rule_a_probe(E: int, stage: int, N: int) -> str:
    """sha256 of the FIRST replicate's resampling indices of a stage: Generator(PCG64(SeedSequence(entropy=E, spawn_key=(stage-1,))))
    .integers(low=0, high=N, size=N, dtype=int64), little-endian int64 bytes."""
    import numpy as np
    child = np.random.SeedSequence(entropy=int(E), spawn_key=(stage - 1,))
    idx = np.random.Generator(np.random.PCG64(child)).integers(low=0, high=N, size=N, dtype=np.int64, endpoint=False)
    return hashlib.sha256(idx.astype("<i8").tobytes()).hexdigest()


def _rule_a_validator(meta, spec, N):
    bad = []
    E = int(spec["root_entropy"])
    if meta["bit_generator"] != spec["bit_generator"] or spec["bit_generator"] != "PCG64":
        bad.append(f"bit generator {meta['bit_generator']!r} is not PCG64")
    if meta["numpy_version"] != spec["numpy_version"]:
        bad.append(f"recorded numpy {meta['numpy_version']!r} is not the pinned {spec['numpy_version']!r}")
    if meta["resampling_api"] != spec["resampling_api"]:
        bad.append("recorded resampling API is not the approved call")
    for stage, rec in ((1, meta["stage1"]), (2, meta["stage2"])):
        if rec is None:
            continue
        sid, seed, probe = rec
        if seed != E:
            bad.append(f"Stage-{stage} seed {seed} is not the recorded root entropy")
        if sid != rule_a_stream_id(E, stage, spec["numpy_version"]):
            bad.append(f"Stage-{stage} stream identifier does not follow the approved format/spawn key")
        if numpy_version() == spec["numpy_version"] and probe != rule_a_probe(E, stage, N):
            bad.append(f"Stage-{stage} stream probe does not reproduce from the recorded metadata")
    return bad


RNG_RULES[RULE_A] = _rule_a_validator


# ----- approved C0 selection (R7): dedicated SeedSequence(entropy=S), PCG64, permutation of the ordered ensemble units -----
C0_ALG = "SEEDSEQ-PCG64-PERMUTATION-C0"


def c0_subset_size(fraction, n_cells, N) -> int:
    """n = ceil(fraction * n_cells * N) in EXACT rational arithmetic (fraction parsed from its decimal text), e.g. (16 N + 99) // 100."""
    return math.ceil(Fraction(str(fraction)) * n_cells * N)


def _c0_select(units, n, seed, spec):
    import numpy as np
    order = sorted(units, key=lambda u: (float(u[0]), float(u[1]), int(u[2])))       # W, zeta, seed ascending; keys unique -> no ties
    perm = np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy=int(seed)))).permutation(len(order))
    return [order[i] for i in sorted(int(i) for i in perm[:n])]


SELECTION_ALGORITHMS[C0_ALG] = _c0_select


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
    # required: Tier 3 at 1/1024 and 1/2048 for the five points; optional (explicitly allowed): the same points at the 1/512 Tier 3
    # validation resolution, where an M4 solve also exists. Any other key is unexpected.
    specs.append(GateSpec("B3Y.decomposition.evaluable", "every M2/M3 decomposition row is evaluable (gates, mirrors, finite terms)",
                          "error_decomposition.csv", ("point", "inv_dx", "method"),
                          expected_keys=[(p, i, m) for p in t3f for i in (mid, fine) for m in ("M2", "M3")],
                          allowed_keys_fn=lambda load: b3y_optional_decomposition_keys(load, t3f, t3v, base),
                          status_column="status", status_ok=("", "OK"), status_not_evaluable=("NOT EVALUABLE",),
                          checks=[Check(c, "finite") for c in ("total", "bulk_rounding", "lead_data", "solver", "identity_residual_three_term", "e_ind_combined")] +
                          [Check(c, "true") for c in ("gate_fl_passed", "gate_hybrid_passed", "mirror_source_bitwise", "mirror_extraction_bitwise",
                                                      "M2_recomputed_bitwise", "M3_recomputed_bitwise")]))
    xg = [("G0", [(p, "") for p in t1]), ("X1", [(p, i) for p in dict.fromkeys(t3v + t3f) for i in inv]), ("X2", [(p, i) for p in t3f for i in inv]),
          ("X3", [(p, i) for p in t3f for i in inv]), ("X6", [(p, i) for p in t3f for i in (base, mid)])]
    for g, keys in xg:
        specs.append(GateSpec(f"B3Y.xp_gate.{g}", f"gate {g}: every protocol-required row passed", "xp_gates.csv", ("case", "inv_dx"), expected_keys=keys,
                              row_filter={"gate": g}, checks=[Check("passed", "true")]))
    for g in ("X0", "X8"):
        specs.append(GateSpec(f"B3Y.xp_gate.{g}", f"gate {g}: exactly one row, passed", "xp_gates.csv", ("gate",), expected_keys=[(g,)],
                              row_filter={"gate": g}, checks=[Check("passed", "true")]))
    specs.append(GateSpec("B3Y.K1_buffer_control", "K1 buffer control rows", "buffer_control.csv", ("point",), expected_keys=[(p,) for p in buf],
                          checks=[Check("K1_passed_dT1", "true"), Check("S_within_geometry_passed", "true"),
                                  Check("abs_dT1_M1", "le", L["dT1_buffer"]), Check("abs_dT1_M3", "le", L["dT1_buffer"])]))
    specs.append(GateSpec("B3Y.H1_solver_spread", "H1 T1 spread rows", "solver_spread.csv", ("point", "inv_dx"), expected_keys=[(p, i) for p in t2 for i in inv],
                          checks=[Check("passed", "true"), Check("T1_spread", "le", L["T1_method_spread_for_solver_independence"])]))
    specs.append(GateSpec("B3Y.O1_rows_present", "O1 convergence rows present (classification only; completeness gate)", "refinement_convergence.csv",
                          ("point", "method"), expected_keys=[(p, m) for p in t2 for m in ("M1", "M2", "M3", "M4")], checks=[Check("observed_order", "finite")]))
    # potential_identity.csv repeats one device record per POINT that uses the device: the multiplicity is protocol-derived
    mult = Counter()
    for p in t1:
        mult[(pts[p]["config"], base, 0.2)] += 1
        mult[(pts[p]["config"], mid, 0.2)] += 1
    for p in t2:
        mult[(pts[p]["config"], fine, 0.2)] += 1
    for p in buf:
        mult[(pts[p]["config"], mid, 0.4)] += 1
    is_p2 = lambda r: "_vs_" in str(r.get("inv_dx", ""))
    specs.append(GateSpec("B3Y.P1_P4_devices", "potential identity and lead guard of every predeclared device record", "potential_identity.csv",
                          ("config", "inv_dx", "buffer_requested"), expected_multiplicity=dict(mult), row_select=lambda r: not is_p2(r),
                          checks=[Check("P1_rgf_device_vs_direct_bitwise", "true"), Check("P4_lead_guard_ok", "true"), Check("V_lead_max", "le", L["lead_tail_guard"])]))
    specs.append(GateSpec("B3Y.P2_cross_resolution", "coincident-slice potential identity across resolutions", "potential_identity.csv",
                          ("config", "inv_dx"), expected_keys=[(c, f"{i}_vs_{base}") for c in sorted({pts[p]["config"] for p in t2}) for i in (mid, fine)],
                          row_select=is_p2, checks=[Check("P2_potential_bitwise", "true"), Check("coincident_x_bitwise", "true")]))
    return specs


def b3y_optional_decomposition_keys(load, t3f, t3v, base):
    """COMPATIBILITY rows (not protocol-required): M2/M3 decomposition at the base resolution 1/512 is accepted only for a point that is in
    BOTH tier3_xp_1024_2048 and tier3_xp_512_validation AND has exactly one tier3_exact_lattice.csv record at inv_dx = 512 with
    gates_passed == True. No other resolution, point or method can enter through this rule."""
    t3 = load("tier3_exact_lattice.csv") or []
    ok = Counter(norm(r.get("point")) for r in t3 if norm(r.get("inv_dx")) == norm(base) and is_true(r.get("gates_passed")))
    return [(p, base, m) for p in t3f if p in t3v and ok.get(norm(p), 0) == 1 for m in ("M2", "M3")]


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


def test_order(D):
    return [(o, c["id"]) for o in D["statistics"]["observables"] for c in D["contrasts"]]


def physical_gate_spec(D: dict, n_seeds: int) -> list:
    S, L = D["seeds"], D["limits"]
    if n_seeds is None or isinstance(n_seeds, bool) or int(n_seeds) != n_seeds or not (S["N_min"] <= int(n_seeds) <= S["N_max"]):
        raise ValueError(f"N = {n_seeds!r} is not a protocol-admissible confirmatory ensemble size [{S['N_min']}, {S['N_max']}]")
    N = int(n_seeds)
    seeds = list(range(S["confirmatory_first"], S["confirmatory_first"] + N))
    E = list(range(D["energies"]["n_energies"]))
    cells = _cells(D, "all")
    C = D["controls"]
    rs = D["row_status_vocabulary"]
    specs = [GateSpec("PHYS.stages.completed", "every predeclared stage completed (no early termination)", "stages.csv", ("stage",),
                      expected_keys=[(s,) for s in D["stages"]], status_column="status", status_ok=("COMPLETED",))]
    specs.append(GateSpec("PHYS.C0.realizations", "C0 per-realization gates: unitarity, reciprocity, channel count, lead tail, V hash",
                          "realizations.csv", ("W", "zeta", "seed", "energy_index"),
                          expected_keys=[(w, z, s, e) for (w, z) in cells for s in seeds for e in E],
                          status_column="status", status_ok=tuple(rs["ok"]), status_not_evaluable=tuple(rs["not_evaluable"]),
                          checks=[Check("unitarity", "le", L["unitarity"]), Check("reciprocity", "le", L["reciprocity"]),
                                  Check("lead_tail_max", "le", L["lead_tail_guard"]), Check("n_open", "eqcol", "n_open_analytic"),
                                  Check("v_hash", "nonempty"), Check("tau", "finite"), Check("T1", "finite")]))
    # C0 reproducibility subset (docs/CAMPAIGN_STATUS_CONTRACT.md section 5): exact rational size rule; dedicated selection seed
    n_sub = c0_subset_size(L["reproducibility_subset_fraction"], len(cells), N)
    ens = [(w, z, s) for (w, z) in cells for s in seeds]
    sel = D.get("reproducibility_selection")
    rep = GateSpec("PHYS.C0.reproducibility", "C0 bitwise-reproducible rerun of the protocol-selected subset: exact units, count and recorded selection",
                   "reproducibility.csv", ("W", "zeta", "seed"), checks=[Check("bitwise_equal", "true"), Check("selection_algorithm", "nonempty"),
                                                                       Check("selection_seed", "nonempty"), Check("selection_numpy_version", "nonempty")])
    alg = None if sel is None else SELECTION_ALGORITHMS.get(sel.get("algorithm"))
    seed = None if sel is None else sel.get("seed")
    boot = D.get("bootstrap_rng") or {}
    if seed is not None:
        if not _seed_ok(str(seed)):
            raise ValueError(f"protocol C0 selection seed {seed!r} is not a non-negative decimal integer")
        if boot.get("root_entropy") is not None and int(seed) == int(boot["root_entropy"]):
            raise ValueError("the C0 selection seed must be dedicated (different from the bootstrap root entropy)")
    if alg is None or seed is None:
        ens_set = {tuple(norm(x) for x in k) for k in ens}

        def subset_rule(load):
            rows = [r for r in (load("reproducibility.csv") or [])]
            keys = [_row_key(r, ("W", "zeta", "seed")) for r in rows]
            inside = [k for k in dict.fromkeys(keys) if k in ens_set]
            return inside if len(inside) == n_sub else inside + [("__protocol_count__", n_sub, len(inside))]     # forces 'missing' if count wrong
        rep.expected_keys_fn = subset_rule
        rep.not_evaluable_reason = ("PROTOCOL DECISION REQUIRED: C0 selection algorithm not specified" if sel is None else
                                    f"no implementation for C0 selection algorithm {sel.get('algorithm')!r}" if alg is None else
                                    "PROTOCOL DECISION REQUIRED: the dedicated C0 selection seed has not been recorded")
    else:
        if numpy_version() != sel.get("numpy_version"):
            rep.expected_keys = [tuple(k) for k in ens][:0]
            rep.expected_keys_fn = lambda load: [_row_key(r, ("W", "zeta", "seed")) for r in (load("reproducibility.csv") or [])]
            rep.not_evaluable_reason = f"installed numpy {numpy_version()} differs from the pinned {sel.get('numpy_version')}: subset not reproducible here"
        else:
            rep.expected_keys = list(alg(ens, n_sub, int(seed), sel))
            if len(set(rep.expected_keys)) != n_sub or not set(rep.expected_keys) <= set(ens):
                raise ValueError("selection algorithm returned an invalid subset")
        rep.row_evaluator = lambda r, sel=sel: None if (str(r.get("selection_algorithm")).strip() == str(sel["algorithm"])
                                                          and _seed_ok(str(r.get("selection_seed"))) and int(r.get("selection_seed")) == int(sel["seed"])
                                                          and str(r.get("selection_numpy_version")).strip() == str(sel["numpy_version"])) \
            else "selection algorithm / seed / numpy version not the protocol's"
    specs.append(rep)
    c2 = C["C2"]
    specs.append(GateSpec("PHYS.C2.solver_cross_check", "C2 M3 vs M1 max|dS| (unadjudicated or adjudicated exceedance is never PASS)",
                          "solver_cross_check.csv", ("W", "zeta", "seed", "energy_index"),
                          expected_keys=[(w, z, s, e) for (w, z) in _cells(D, c2["cells"]) for s in seeds[:c2["n_seeds"]] for e in E],
                          checks=[Check("max_abs_dS", "le", L["solver_cross_check_S_abs"])]))
    lim = L["solver_cross_check_S_abs"]

    def flagged(load):
        rows = load("solver_cross_check.csv") or []
        out = []
        for r in rows:
            x = num(r.get("max_abs_dS"))
            if x is None or not (math.isfinite(x) and x <= lim):
                out.append(_row_key(r, ("W", "zeta", "seed", "energy_index")))
        return list(dict.fromkeys(out))
    # C2 adjudication, two separate questions (docs/CAMPAIGN_STATUS_CONTRACT.md, section 4):
    #  (a) administrative: exactly one adjudication record per flagged C2 row, with a valid status and the recorded M4 comparison numbers;
    #  (b) scientific: is the discrepancy RESOLVED? Requires a protocol-approved numerical criterion. None is approved, so (b) is NOT
    #      EVALUABLE whenever anything is flagged. Neither gate can turn a flagged C2 row into a PASS: PHYS.C2.solver_cross_check stays FAIL.
    adj_num = [Check("eps_M1", "finite"), Check("eps_M3", "finite")]
    specs.append(GateSpec("PHYS.C2.adjudication_recorded", "administrative: exactly one adjudication record per flagged C2 row (status, M4 numbers)",
                          "adjudication.csv", ("W", "zeta", "seed", "energy_index"), expected_keys_fn=flagged, status_column="status",
                          status_ok=("ADJUDICATED",), status_not_evaluable=("NOT EVALUABLE",), checks=adj_num))
    crit = D.get("C2_adjudication_criteria")
    crit_fn = None if crit is None else C2_ADJUDICATION_CRITERIA.get(crit.get("criterion"))

    def resolution_open(load):
        if not flagged(load):
            return None                                                      # nothing to resolve
        if crit is None:
            return "PROTOCOL DECISION REQUIRED: no numerical C2 adjudication (resolution) criterion is approved"
        if crit_fn is None:
            return f"no implementation for the C2 adjudication criterion {crit.get('criterion')!r}"
        return None
    specs.append(GateSpec("PHYS.C2.adjudication_resolved", "scientific: every flagged C2 discrepancy satisfies the approved numerical resolution criterion",
                          "adjudication.csv", ("W", "zeta", "seed", "energy_index"), expected_keys_fn=flagged, status_column="status",
                          status_ok=("ADJUDICATED",), status_not_evaluable=("NOT EVALUABLE",), checks=adj_num,
                          row_evaluator=(lambda r: crit_fn(r, crit)) if crit_fn is not None else None, not_evaluable_fn=resolution_open))
    c5 = C["C5"]
    specs.append(GateSpec("PHYS.C5.b3y_reproduction", "C5 seed 0 reproduces B3Y", "b3y_reproduction.csv", ("zeta", "energy_index"),
                          expected_keys=[(z, e) for z in c5["zeta"] for e in E], checks=[Check("abs_dT1", "le", L["b3y_reproduction_abs_dT1"])]))
    specs.append(GateSpec("PHYS.C6.clean_guide", "C6 clean guide T_nn = 1", "clean_guide.csv", ("energy_index",), expected_keys=[(e,) for e in E],
                          checks=[Check("max_abs_Tnn_minus_1", "le", L["clean_guide_Tnn_abs"])]))
    tost = tost_evaluator(L["equivalence_fraction_of_h"], L["tost_alpha"], L["underpowered_fraction_of_eps"], L["control_seed_extension_max"])
    contrasts = D["contrasts"]
    for cid in ("C1", "C1b", "C3", "C4"):
        c = C[cid]
        cs_ = set(_cells(D, c["cells"]))
        keys = []
        for o in c["observables"]:
            keys += [(cid, o, "cell", cell_id(w, z)) for (w, z) in sorted(cs_)]
            keys += [(cid, o, "contrast", k["id"]) for k in contrasts if all(tuple(x) in cs_ for x in k["cells"])]
        specs.append(GateSpec(f"PHYS.{cid}.equivalence", f"{cid} {c['what']}: TOST equivalence at alpha = {L['tost_alpha']} with power check",
                              "controls.csv", ("control_id", "observable", "unit_type", "unit_id"), expected_keys=keys,
                              row_filter={"control_id": cid}, row_evaluator=tost))
    c7 = C["C7"]
    keys = [(c7["observable"], "cell", cell_id(w, z)) for (w, z) in cells] + [(c7["observable"], "contrast", k["id"]) for k in contrasts]
    specs.append(GateSpec("PHYS.C7.precision", "C7 achieved 95 % CI half-width <= delta (primary observable only)", "precision.csv",
                          ("observable", "unit_type", "unit_id"), expected_keys=keys,
                          checks=[Check("half_width", "le", L["delta"])]))
    order = test_order(D)
    if len(order) != D["statistics"]["n_tests"]:
        raise ValueError("draft: number of tests inconsistent with the contrast list")
    specs.append(GateSpec("PHYS.statistics.tests_complete", "68 confirmatory tests: section 5.3 schema, MC intervals, Stage-2 selection, Holm sets, decisions",
                          "tests.csv", ("observable", "contrast_id"), expected_keys=order,
                          status_column="decision", status_ok=tuple(D["statistics"]["decisions"]),
                          table_evaluator=statistics_table_evaluator(D, order)))
    specs.append(bootstrap_rng_gate(D, order, N))
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
        print(f"{g['state']:14s} {g['gate_id']:44s} expected={g['n_expected']} rows={g['n_rows']} {g['reason']}")
    print("OVERALL:", r["overall_status"], r["problems"])
