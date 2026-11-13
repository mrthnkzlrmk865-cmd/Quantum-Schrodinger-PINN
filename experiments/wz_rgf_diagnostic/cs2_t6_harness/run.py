"""Staged T1-T6 software-validation runner (PREDICTIONS_coupling_scaling_v2.md §2.4.3.11 order of execution).

Stages, in order (a later stage runs only if every earlier stage passed; the official run stops at the first failure):
  1. environment and import checks (mpmath == 1.3.0; AST import guard, self-tested);
  2. clean-room reference known-answer tests and self-proofs (cs2_t6_reference), including the P5 N1' self-check;
  3. frozen-input preconditions (P5 N1-N5 and status list, mutation witnesses, targeted-lattice sensitivity, P4 omega);
  4. T1-T5 (§2.4.2; package cs2_t1t5, against coupling_scaling_v2_certify and coupling_scaling_v2_covering);
  5. T6.1-T6.8 (exact), then the C4 field "numerical cross-check (not validated)" (only its blocking outcomes gate) and
     the diagnostic double-precision anchors (never part of PASS);
  6. the mutation suite MU1-MU27 (36 variants including the supplementary MU18-N4).

Exit codes: 0 = official run, every stage passed; 1 = a stage failed (the run stopped there); 3 = development run.
`--dev-skip-t1t5` is a development aid: stages 5-6 run even if stage 4 fails, the run is labelled DEVELOPMENT, it never
exits 0 and it is never an official result. `--no-c4` skips the C4 field and records a report note. CAUTION (flagged for R0): the exit-code
logic does not distinguish such a run, so `--no-c4` with all other stages passing still exits 0; a run without C4 is
NOT a complete official validation and must not be cited as one (the official command line has no flags).
This runner does NOT generate the manifest and is NOT the R0 review (§8.2).
Usage (from experiments/wz_rgf_diagnostic):  python -m cs2_t6_harness.run --out <dir> [--dev-skip-t1t5] [--no-c4]
"""
import argparse
import sys
import time

import mpmath

import coupling_scaling_v2_certify as prod
from .common import Results
from . import import_guard, kat, witnesses, t6_1, t6_2, t6_3, t6_4, t6_5, t6_6, t6_7, t6_8, pointwise_cache, c4, anchors
from . import mutations, report
import coupling_scaling_v2_covering as cov
from cs2_t1t5 import run_all as t1t5


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--dev-skip-t1t5", action="store_true")
    ap.add_argument("--no-c4", action="store_true")
    a = ap.parse_args(argv)
    stages, notes = [], []
    payload = dict(mode="DEVELOPMENT (stage-4 result ignored; not an official T6 result)" if a.dev_skip_t1t5 else "OFFICIAL",
                   stages=stages, notes=notes)

    def record(R, status=None, t0=None):
        st = R.summary()
        st["status"] = status or ("PASS" if not R.failures else "FAIL")
        st["seconds"] = round(time.time() - t0, 1) if t0 else None
        st["failures"] = R.failures[:50]
        stages.append(st)
        print(f"[stage] {st['stage']}: {st['total']} checks, {st['failed']} failed -> {st['status']}", flush=True)
        return not R.failures

    def finish(overall):
        payload["overall"] = overall
        out = report.write(a.out, payload)
        print(f"[overall] {overall}  (report: {out})")
        # Exit 0 only for an official run (no --dev-skip-t1t5) in which every stage 1-6 passed.
        # A development run never exits 0: 3 = all development stages passed (stage-4 result not gating); 1 = a stage
        # failed and the run stopped there.
        if a.dev_skip_t1t5:
            return 3 if overall.startswith("ALL STAGES PASS") else 1
        return 0 if overall.startswith("ALL STAGES PASS") else 1

    t0 = time.time()
    R = Results("1 environment and import checks")
    R.add(f"mpmath == 1.3.0 (found {mpmath.__version__})", mpmath.__version__ == "1.3.0", "ENV-MISMATCH")
    import_guard.check(R)
    if not record(R, t0=t0):
        return finish("STOPPED at stage 1 (ENV-MISMATCH / IMPORT-VIOLATION)")
    t0 = time.time()
    if not record(kat.run(), t0=t0):
        return finish("STOPPED at stage 2 (KAT-FAIL / REF-SELFCHECK-FAIL)")
    t0 = time.time()
    if not record(witnesses.run(), t0=t0):
        return finish("STOPPED at stage 3 (PRECONDITION-UNMET)")
    t0 = time.time()
    R = Results("4 T1-T5 (§2.4.2)")
    sub = []
    for Rt in t1t5.run(prod, cov):
        sub.append(Rt.summary())
        for it in Rt.items:
            R.items.append(it)
    payload["t1t5"] = sub
    for st in sub:
        print(f"        {st['stage']}: {st['total']} checks, {st['failed']} failed", flush=True)
    ok4 = record(R, t0=t0)
    if not ok4 and not a.dev_skip_t1t5:
        return finish("STOPPED at stage 4 (T1-T5 failure)")
    if a.dev_skip_t1t5:
        notes.append("DEVELOPMENT run (--dev-skip-t1t5): stages 5-6 run regardless of stage 4. Not an official T6 result.")
    t0 = time.time()
    R = Results("5 T6.1-T6.8")
    t6_1.run(prod, R)
    t6_2.run(prod, R)
    t6_3.run(prod, R)
    t6_4.run(prod, R)
    cache = pointwise_cache.compute(R, prod)
    t6_5.run(prod, R, cache=cache)
    t6_6.run(prod, R, cache=cache)
    t6_7.run(prod, R)
    t6_8.run(prod, R)
    ok5 = record(R, t0=t0)
    if not a.no_c4:
        t0 = time.time()
        Rc, info = c4.run(prod)
        payload["c4"] = info
        st = Rc.summary()
        st["status"] = ("BLOCKING " if info["blocking"] else "non-blocking ") + info["field_outcome"]
        st["seconds"] = round(time.time() - t0, 1)
        st["failures"] = Rc.failures[:50]
        stages.append(st)
        print(f"[stage] C4 numerical cross-check (not validated): field outcome {info['field_outcome']}, "
              f"blocking={info['blocking']}", flush=True)
        ok5 = ok5 and not info["blocking"]
    else:
        notes.append("C4 not run (--no-c4).")
    payload["anchors"] = anchors.run()
    notes.append("Double-precision anchors are diagnostic only and never part of PASS.")
    if not ok5:
        return finish("STOPPED at stage 5 (T6 failure or blocking C4 outcome)")
    t0 = time.time()
    Rm, rows = mutations.run()
    payload["mutations"] = rows
    if not record(Rm, t0=t0):
        return finish("STOPPED at stage 6 (MUTATION-SURVIVED)")
    if a.dev_skip_t1t5:
        return finish("ALL STAGES PASS IN DEVELOPMENT MODE (stage-4 result not gating) — NOT AN OFFICIAL T6 PASS")
    return finish("ALL STAGES PASS (OFFICIAL): T1-T6 software self-tests; manifest NOT generated by this runner")


if __name__ == "__main__":
    sys.exit(main())
