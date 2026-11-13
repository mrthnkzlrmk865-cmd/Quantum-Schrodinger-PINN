"""Fast T6 harness checks for the repository test suite: import guard, reference KATs/self-proofs, frozen-input
preconditions (incl. P5 N1-N5), and T6.7. The full staged run is `python -m cs2_t6_harness.run` (see that module)."""
import sys
import pathlib

import pytest

HERE = pathlib.Path(__file__).resolve().parent.parent / "experiments" / "wz_rgf_diagnostic"
sys.path.insert(0, str(HERE))

mpmath = pytest.importorskip("mpmath")


def test_import_guard():
    from cs2_t6_harness import import_guard
    from cs2_t6_harness.common import Results
    R = import_guard.check(Results("guard"))
    assert not R.failures, R.failures


def test_reference_kats_and_selfproofs():
    from cs2_t6_harness import kat
    R = kat.run()
    assert not R.failures, R.failures[:5]


def test_frozen_preconditions_including_p5():
    from cs2_t6_harness import witnesses
    R = witnesses.run()
    assert not R.failures, R.failures[:5]


def test_t6_7_threshold_and_comparison():
    import coupling_scaling_v2_certify as prod
    from cs2_t6_harness import t6_7
    R = t6_7.run(prod)
    assert not R.failures, R.failures


def test_t1_to_t5():
    import coupling_scaling_v2_certify as prod
    import coupling_scaling_v2_covering as cov
    from cs2_t1t5 import run_all
    results = run_all.run(prod, cov)
    assert len(results) == 5
    for R in results:
        assert not R.failures, (R.stage, R.failures[:5])
