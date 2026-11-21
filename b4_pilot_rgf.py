"""B4 limited diagnostic pilot — RGF order classification and deterministic resolution selection (phase T3).

Classification: DIAGNOSTIC / NOT B4.

Pure functions only: no file I/O, no git, no solver execution, no result-directory writes. They implement
docs/B4_PILOT_PROTOCOL.md §1.4 (observed-order status) and §1.6 (admissibility and deterministic selection of the
RGF resolution) exactly as written there. No threshold is defined here that the protocol does not already state.

A SELECTED resolution is an input to gate G-RGF, not a G-RGF pass: G-RGF also needs solver agreement on physical
data, the N_y/N_q checks and the lead-tail guard (protocol §6). No physical RGF data are evaluated in this module.

Richardson values are estimates of the dx -> 0 limit of the finite-difference problem at fixed N_y and N_q under the
asymptotic assumption. They are not demonstrated continuum convergence. BELOW FLOOR / ADMISSIBLE-F is not an error
bound.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping, Sequence

CLASSIFICATION = "DIAGNOSTIC / NOT B4"

# ---- protocol values (§0.1, §1.4, §1.6); quoted, never tuned ----------------------------------------------------
ORDER_FLOOR = 1e-8                 # §0.1 / §1.4: a change is "below the floor" iff |Δ| < 1e-8 (strict)
ORDER_WINDOW = (1.8, 2.2)          # §1.4: RESOLVED iff p in [1.8, 2.2] (closed)
CHANGE_LIMIT = 2e-4                # §1.6: consecutive-change limit, |Δ| <= 2e-4 (inclusive)
ERROR_LIMIT = 1e-4                 # §1.6: Richardson-estimate limit, e <= 1e-4 (inclusive)

DX_CANDIDATES = ("1/512", "1/1024", "1/2048")                 # coarse -> fine
RICHARDSON_K = {"1/512": 2, "1/1024": 1, "1/2048": 0}          # §1.6: e(dx_c) = 2^{k p} |Δ₂| / (2^p - 1)
QUANTITIES = ("T1", "R1")

# §1.4 statuses
BELOW_FLOOR = "BELOW FLOOR"
NOT_RESOLVABLE = "NOT RESOLVABLE"
OUTSIDE_WINDOW = "OUTSIDE WINDOW"
RESOLVED = "RESOLVED"

# §1.6 labels. INADMISSIBLE (a RESOLVED quantity that fails a limit) is this module's name for the protocol's
# "not ADMISSIBLE"; the protocol table names only the four others.
ADMISSIBLE = "ADMISSIBLE"
ADMISSIBLE_R = "ADMISSIBLE-R"      # dx = 1/2048 only: Richardson-only, no next-finer change available
ADMISSIBLE_F = "ADMISSIBLE-F"      # below floor: NOT an error bound
INADMISSIBLE_U = "INADMISSIBLE-U"  # unclassifiable (NOT RESOLVABLE or OUTSIDE WINDOW)
INADMISSIBLE = "INADMISSIBLE"
ADMISSIBLE_LABELS = frozenset({ADMISSIBLE, ADMISSIBLE_R, ADMISSIBLE_F})

LABEL_MEANING = {
    ADMISSIBLE: "RESOLVED; change to the next finer grid <= 2e-4 and Richardson estimate <= 1e-4",
    ADMISSIBLE_R: "RESOLVED at dx = 1/2048: no next-finer change exists (no 1/4096); admissible on the change into "
                  "1/2048 <= 2e-4 and the Richardson estimate <= 1e-4 under the observed-order assumption",
    ADMISSIBLE_F: "BELOW FLOOR: both changes < 1e-8; not an error bound, not an asymptotic regime, not continuum "
                  "convergence",
    INADMISSIBLE_U: "NOT RESOLVABLE or OUTSIDE WINDOW: the finest available pair cannot classify the order",
    INADMISSIBLE: "RESOLVED but a limit (2e-4 change or 1e-4 Richardson estimate) is exceeded",
}

# selection outcomes
SELECTED = "SELECTED"
NOT_EVALUABLE = "G-RGF NOT EVALUABLE"
FAIL = "G-RGF FAIL"
STOP = "STOP: SOLVER AGREEMENT NOT ESTABLISHED"


@dataclass(frozen=True)
class OrderResult:
    """§1.4 result for one quantity at one point. p is None unless both changes are >= floor with the same sign."""
    q512: float
    q1024: float
    q2048: float
    delta1: float                  # Q(1/1024) - Q(1/512)
    delta2: float                  # Q(1/2048) - Q(1/1024)
    status: str
    p: float | None
    classification: str = CLASSIFICATION


@dataclass(frozen=True)
class CandidateResult:
    dx: str
    label: str
    error_estimate: float | None   # e(dx) for RESOLVED; None otherwise (never fabricated)
    change: float | None           # |Δ| used by the change condition for RESOLVED; None otherwise


@dataclass(frozen=True)
class SelectionResult:
    outcome: str                                      # SELECTED, NOT_EVALUABLE, FAIL or STOP
    dx_selected: str | None
    orders: dict = field(default_factory=dict)        # (point, quantity) -> OrderResult
    candidates: dict = field(default_factory=dict)    # (point, quantity) -> {dx: CandidateResult}
    labels_at_selected: dict = field(default_factory=dict)   # (point, quantity) -> label at dx_selected
    status_mismatch_points: tuple = ()                # points where T1 and R1 statuses differ (§1.4: flag)
    statements: tuple = ()                            # §1.6 step 4 explicit statements
    classification: str = CLASSIFICATION


def _finite(x, name):
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        raise TypeError(f"{name} must be a real number, got {type(x).__name__}")
    x = float(x)
    if not math.isfinite(x):
        raise ValueError(f"{name} must be finite, got {x!r}")
    return x


def classify_order(q512: float, q1024: float, q2048: float) -> OrderResult:
    """§1.4 status of ONE quantity from its values at dx = 1/512, 1/1024, 1/2048.

    Δ₁ = Q(1/1024) - Q(1/512), Δ₂ = Q(1/2048) - Q(1/1024).
    BELOW FLOOR     |Δ₁| < 1e-8 and |Δ₂| < 1e-8
    NOT RESOLVABLE  exactly one of |Δ₁|, |Δ₂| < 1e-8, or sign(Δ₁) != sign(Δ₂)
    OUTSIDE WINDOW  both >= 1e-8, same sign, p = log2(|Δ₁|/|Δ₂|) not in [1.8, 2.2]
    RESOLVED        both >= 1e-8, same sign, p in [1.8, 2.2] (endpoints included)
    p is returned only when it is defined (both >= floor, same sign); otherwise None.
    """
    q512, q1024, q2048 = (_finite(q512, "q512"), _finite(q1024, "q1024"), _finite(q2048, "q2048"))
    d1 = q1024 - q512
    d2 = q2048 - q1024
    below1, below2 = abs(d1) < ORDER_FLOOR, abs(d2) < ORDER_FLOOR
    if below1 and below2:
        return OrderResult(q512, q1024, q2048, d1, d2, BELOW_FLOOR, None)
    if below1 or below2:
        return OrderResult(q512, q1024, q2048, d1, d2, NOT_RESOLVABLE, None)
    if (d1 > 0) != (d2 > 0):                                   # both nonzero here (>= floor)
        return OrderResult(q512, q1024, q2048, d1, d2, NOT_RESOLVABLE, None)
    p = math.log2(abs(d1) / abs(d2))
    lo, hi = ORDER_WINDOW
    status = RESOLVED if lo <= p <= hi else OUTSIDE_WINDOW
    return OrderResult(q512, q1024, q2048, d1, d2, status, p)


def richardson_error(order: OrderResult, dx: str) -> float:
    """§1.6 e(dx_c) = 2^{k p} |Δ₂| / (2^p - 1), k = 2, 1, 0 for dx_c = 1/512, 1/1024, 1/2048.
    Defined only for a RESOLVED quantity; raises ValueError otherwise (no estimate is fabricated)."""
    if dx not in RICHARDSON_K:
        raise ValueError(f"unknown candidate dx {dx!r}; expected one of {DX_CANDIDATES}")
    if order.status != RESOLVED or order.p is None:
        raise ValueError(f"Richardson estimate is NOT APPLICABLE for status {order.status!r}")
    p = order.p
    return 2.0 ** (RICHARDSON_K[dx] * p) * abs(order.delta2) / (2.0 ** p - 1.0)


def admissibility(order: OrderResult) -> dict:
    """§1.6 step 2: label of every candidate dx for one (point, quantity). Returns {dx: CandidateResult}."""
    out = {}
    for dx in DX_CANDIDATES:
        if order.status == BELOW_FLOOR:
            out[dx] = CandidateResult(dx, ADMISSIBLE_F, None, None)
        elif order.status in (NOT_RESOLVABLE, OUTSIDE_WINDOW):
            out[dx] = CandidateResult(dx, INADMISSIBLE_U, None, None)
        elif order.status == RESOLVED:
            e = richardson_error(order, dx)
            change = abs(order.delta1) if dx == "1/512" else abs(order.delta2)
            ok = change <= CHANGE_LIMIT and e <= ERROR_LIMIT
            good = ADMISSIBLE_R if dx == "1/2048" else ADMISSIBLE
            out[dx] = CandidateResult(dx, good if ok else INADMISSIBLE, e, change)
        else:                                                  # pragma: no cover - defensive
            raise ValueError(f"unknown status {order.status!r}")
    return out


def select_resolution(values: Mapping, required_points: Sequence, *, solver_agreement: bool) -> SelectionResult:
    """§1.6 steps 0-4 for ONE configuration.

    values           mapping (point, quantity) -> (Q(1/512), Q(1/1024), Q(1/2048)), quantity in ("T1", "R1"),
                     M1 values. Every required point must have BOTH quantities; missing, extra or malformed entries
                     raise ValueError (nothing is silently omitted; R1 is never inferred from T1).
    required_points  the complete §1.6 point set of the configuration (no duplicates, non-empty).
    solver_agreement step 0: True only if M1-M2-M3 agreement (<= 1e-10, or adjudicated as in B3Y) holds at every
                     point and dx. Mandatory keyword; False -> outcome STOP.

    Outcome: SELECTED with dx_selected = coarsest candidate admissible (ADMISSIBLE, ADMISSIBLE-R or ADMISSIBLE-F)
    for every (point, quantity); otherwise NOT EVALUABLE if any (point, quantity) is INADMISSIBLE-U, else FAIL.
    """
    if not isinstance(solver_agreement, bool):
        raise TypeError("solver_agreement must be given explicitly as True or False")
    points = list(required_points)
    if not points:
        raise ValueError("required_points is empty")
    if len(set(points)) != len(points):
        raise ValueError("required_points contains duplicates")
    expected = {(x, q) for x in points for q in QUANTITIES}
    given = set(values.keys())
    missing = sorted(expected - given, key=repr)
    extra = sorted(given - expected, key=repr)
    if missing or extra:
        raise ValueError(f"incomplete or inconsistent inputs: missing {missing}; unexpected {extra}")
    if not solver_agreement:
        return SelectionResult(STOP, None)

    orders, cands = {}, {}
    for x in points:
        for q in QUANTITIES:
            v = values[(x, q)]
            if v is None or len(v) != 3:
                raise ValueError(f"{(x, q)}: expected three values (dx = 1/512, 1/1024, 1/2048), got {v!r}")
            orders[(x, q)] = classify_order(*v)
            cands[(x, q)] = admissibility(orders[(x, q)])

    mismatch = tuple(x for x in points if orders[(x, "T1")].status != orders[(x, "R1")].status)
    keys = [(x, q) for x in points for q in QUANTITIES]
    for dx in DX_CANDIDATES:                                   # coarse -> fine: first all-admissible is the coarsest
        if all(cands[k][dx].label in ADMISSIBLE_LABELS for k in keys):
            labels = {k: cands[k][dx].label for k in keys}
            st = []
            if dx == "1/2048":
                st.append("dx_sel = 1/2048: no next-finer change is available; admissibility rests on the "
                          "observed-order assumption from three grids")
            if any(lab == ADMISSIBLE_F for lab in labels.values()):
                st.append("at least one (point, quantity) is ADMISSIBLE-F: below floor, NOT an error bound")
            return SelectionResult(SELECTED, dx, orders, cands, labels, mismatch, tuple(st))
    if any(cands[k][dx].label == INADMISSIBLE_U for k in keys for dx in DX_CANDIDATES):
        return SelectionResult(NOT_EVALUABLE, None, orders, cands, {}, mismatch,
                               ("the finest available pair cannot classify the order for at least one "
                                "(point, quantity); remedy: dated pilot amendment adding dx = 1/4096 (§1.5)",))
    return SelectionResult(FAIL, None, orders, cands, {}, mismatch,
                           ("every (point, quantity) is classified but a limit is exceeded even at dx = 1/2048",))
