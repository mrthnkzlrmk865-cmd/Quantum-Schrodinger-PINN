import os
import sys

# Performance setting for the test session (must precede the first NumPy
# import). With a multithreaded OpenBLAS, the small (63x63) complex BLAS
# calls of the RGF solver suffered ~4-16 ms thread wake-up latency each on
# the 4-core validation machine (~30-40x slower). For the clean-guide RGF
# results, 1, 2 and 4 threads gave bit-identical S matrices
# (docs/RGF_B1_VALIDATION.md). setdefault: an explicit user setting wins.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
