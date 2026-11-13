"""Clean-room T6 reference (PREDICTIONS_coupling_scaling_v2.md §2.4.3.4).

Exact rational arithmetic only (fractions.Fraction, math.isqrt). This package must not import the production certificate
module, the harness, wz_mutants, coupled_blocks, wz_rgf, stationary, numpy, or any other repository module. Only `rnum`
(the C4 numerical cross-check) may import mpmath. The rules are enforced by cs2_t6_harness/import_guard.py.
"""
