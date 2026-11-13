"""Staged T6 runner (§2.4.3.11 order of execution):
  1 ENV and import checks; 2 reference KATs and self-proofs; 3 frozen-input preconditions; 4 T1-T5; 5 T6.1-T6.8 (+ C4 field,
  diagnostic anchors); 6 mutation suite. A later stage runs only if every earlier stage passed.

T1-T5 are NOT implemented (§2.4.2 / §8.1), so stage 4 is BLOCKED and the official run stops there. `--dev-skip-t1t5`
continues past stage 4 for development; that run is labelled DEVELOPMENT and is never an official T6 result.
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--dev-skip-t1t5", action="store_true")
    ap.add_argument("--no-c4", action="store_true")
    a = ap.parse_args(argv)
    stages, notes = [], []
    payload = dict(mode="DEVELOPMENT (stage 4 bypassed; not an official T6 result)" if a.dev_skip_t1t5 else "OFFICIAL",
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
        # Exit 0 only for an official run with every stage passing (impossible while T1-T5 are not implemented).
        # A development run never exits 0: 3 = development stages passed, stage 4 BLOCKED; 1 = a stage failed/stopped.
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
    R = Results("4 T1-T5 (§2.4.2)")
    R.add("T1-T5 implemented", False, "BLOCKED", "T1-T5 are specified in §2.4.2 but not implemented (§8.1 planned files)")
    record(R, status="BLOCKED")
    notes.append("Stage 4 (T1-T5) is BLOCKED: not implemented. The manifest-freeze condition 'T1-T5 pass' cannot be met.")
    if not a.dev_skip_t1t5:
        return finish("STOPPED at stage 4 (T1-T5 BLOCKED: not implemented)")
    notes.append("DEVELOPMENT run: stages 5-6 executed despite stage 4 being BLOCKED. Not an official T6 result.")
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
    return finish("ALL STAGES PASS EXCEPT STAGE 4 (BLOCKED) — DEVELOPMENT RUN, NOT AN OFFICIAL T6 PASS")


if __name__ == "__main__":
    sys.exit(main())
