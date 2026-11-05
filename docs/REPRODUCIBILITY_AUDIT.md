# Phase 1 audit: reproducibility, cache correctness, baseline verification

Date: 2026-11-05. Starting point: `origin/main` at `4e24e9e` (merge of PR #1).
Scope: make the existing research computationally reproducible. **No change
to the physical model, the Hamiltonian, the parameters, the PINN
architecture or any scientific conclusion.** No new parameter sweeps.

## 1. Repository protection

* The branch `claude/quantum-schrodinger-audit-iy6kn4` had no uncommitted
  changes. `origin/main` was this branch plus the PR #1 merge commit, with
  identical file contents. The branch was **fast-forwarded** onto
  `origin/main` (no reset, rebase or force push).
* All committed results under `results/` (validation, transport,
  localization, PINN) and all source files are unchanged in content. New
  results were written only to the fresh directory `results/reproduction/`.
* Physics source files (`config.py`, `disorder.py`, `potential.py`,
  `initial_condition.py`, `reference_solver.py`, `observables.py`,
  `transport.py`) differ from `origin/main` by **0 lines**.

## 2. Cache defects found

| Driver | Reuse rule before | Defect |
|---|---|---|
| `experiment_disorder.py` | file exists and has ≥ n_real realizations | ignored every configuration parameter (W, ζ, grid, Δt, t_end, …); silently returned *more* realizations than requested |
| `experiment_localization.py` | file exists | ignored configuration **and** the realization count |
| `experiment_pinn.py` | `model.pt` exists | ignored configuration; `--quick` runs used the same directory as full runs, so a full run could silently reuse a quick model |
| all | – | no way to write to a fresh directory; a rerun would have read (or partly overwritten) the committed baseline |

Metadata defect: in **22 of 44** committed transport `.npz` files, the stored
`config_json` says `T = 1.0`, while the propagation actually ran to
t_end = 1.25–2.5 (barrier, length and slow-momentum points). The
scientific numbers are **not affected**: `summary.csv` and all analysis
take t_end from the stored `times` array. But the old metadata alone could
not reproduce those runs. The new metadata stores `run.t_end` explicitly.

## 3. Changes (files and reasons)

| File | Change | Reason |
|---|---|---|
| `provenance.py` (new) | Canonical configuration **fingerprint** (SHA-256 of sorted JSON; numbers compared by value), **code fingerprint** (SHA-256 of the source files that determine each experiment's numbers, plus per-file hashes), run metadata (spec, seeds, n_real, UTC time, git commit and dirty flag, Python/NumPy/SciPy/PyTorch versions, platform, CPU count, wall time), and `check_cache()` | single, tested cache policy |
| `experiment_disorder.py` | fingerprinted cache; spec = solver-relevant config + t_end + n_frames + explicit seed list + n_real; metadata stored inside each `.npz` (`metadata_json`, `fingerprint`); sweep-level `run_metadata.json`; `fingerprint` column in `summary.csv`; `--out`, `--recompute`, `--strict-code` | cache correctness, provenance |
| `experiment_localization.py` | same (`run_or_load`), seed list includes the deterministic single W = 0 run | same |
| `experiment_pinn.py` | fingerprint of the full PINN config (excluding `device` and `log_every`); `run_metadata.json` per run (records torch thread count); `--quick` runs go to `<run>_quick`; `--out`, `--recompute`, `--strict-code` | same |
| `validate_solver.py` | `--out`; writes `run_metadata.json` | provenance; fresh directory |
| `compare_reproduction.py` (new) | array-by-array comparison of a rerun with the committed baseline, plus trend checks | baseline verification |
| `tests/test_cache.py` (new) | 39 tests (Section 4) | required tests |
| `requirements-lock.txt` (new) | exact library versions used for the committed and reproduced results | bit-level reproducibility |
| `README.md` | reproduction section: cache policy, new options, lock file, link to this report | documentation |
| `results/reproduction/` (new) | reproduction outputs with full metadata, `resources.jsonl` (wall time, CPU time, peak memory), `comparison.json`, `run_reproduction.sh` | baseline verification |

**Cache policy.** A cached result is reused only if its configuration
fingerprint equals that of the request. A mismatching cache, or a legacy
cache without metadata, raises `CacheMismatchError`. It is never used
silently and never overwritten unless `--recompute` is given. If only the
source code changed, a warning is printed (`--strict-code` makes it an
error). Fields that do not change the numbers are excluded from the
fingerprint and recorded as metadata: device, log interval, batch size,
FFT workers, torch threads. Note that the PINN *training* result does depend
on the torch thread count (Section 5); the count is recorded, but it is not
part of the fingerprint.

## 4. Tests executed (actual outcomes)

Command (from the README): `python -m pytest -q tests`

| Run | Result |
|---|---|
| untouched `origin/main` (temporary worktree) | **22 passed, 1 xfailed** in 13 s |
| this branch | **61 passed, 1 xfailed** in 28 s (22 existing + 39 new) |

Pre-existing outcome (unchanged): the one `xfail` is
`test_trained_pinn_meets_20_percent_accuracy`, which documents that the
PINN does not reach 20 % error (README §9.10). There were no pre-existing
failures, and the changes introduced no regressions.

New cache tests (`tests/test_cache.py`), all passing:
* fingerprint canonicalization (key order, tuple/list, 100 vs 100.0) and
  sensitivity to a 1-ulp change;
* transport spec sensitive to each of 13 physical/solver fields, to the
  realization count, t_end and n_frames, and *insensitive* to PINN-only
  fields; PINN spec sensitive to 9 training/model/problem fields and
  insensitive to device and logging;
* real (tiny) transport runs: reuse without recomputation and with
  identical arrays; rejection after a change of W, ζ, n_real (more *or*
  fewer) or t_end, with the cached file left byte-identical;
  `--recompute` overwrites; legacy files without metadata are rejected;
  a code change warns, or fails with `strict_code`; metadata completeness;
  sweep summary and directory metadata;
* real (tiny) localization runs: reuse, rejection on n_real and t_end
  changes;
* PINN runs (training stubbed): reuse, rejection on config change, quick
  runs isolated, legacy run directories rejected.

End-to-end check on real data: re-invoking the full strength sweep on
`results/reproduction/transport` gave 8/8 cache hits in 2.9 s (vs 15 min),
no code warning, and a byte-identical `summary.csv`. Requesting
`--n-real 40` was refused. Running with the default output directory
(the committed legacy baseline) was refused as unverifiable, and
`git status` shows that nothing in `results/transport` was modified.

## 5. Baseline reproduction

Script: `results/reproduction/run_reproduction.sh` (plus one extra PINN
run with 2 threads). Comparison: `python compare_reproduction.py` →
`results/reproduction/comparison.json`. Large sweeps were not rerun.

| Experiment | Size | Result vs committed baseline |
|---|---|---|
| Solver validation V1–V7 | full | **bit-identical** (all 6 CSVs, full float precision) |
| Transport, strength sweep | 8 points × 48 realizations | **bit-identical**: all raw per-realization T, R, P_in arrays and `summary.csv` (max \|diff\| = 0) |
| Transport, correlation-length sweep | 7 points × 48 | **bit-identical** (max \|diff\| = 0) |
| Closed-box localization | 8 points × 16 | **bit-identical** (MSD, PR arrays and summary) |
| PINN `clean_fourier_uniform_s0`, 2 torch threads (as in the baseline) | 1 run | trained **weights bit-identical**; evaluation metrics equal to ≤ 7.5×10⁻⁷ (see below) |
| PINN `clean_fourier_uniform_s0`, 4 torch threads | 1 run | weights differ (max \|Δw\| = 0.011); final rel. L2 0.48850 → 0.48987, mean 0.31300 → 0.31357, \|ψ\|² 0.5736 → 0.5811 |

Trends: in both baseline and reproduction, ⟨T⟩ decreases monotonically
with W, and the minimum of ⟨T⟩ over ζ is at ζ = 0.05 (kζ ≈ 1). These, and
every number in README §9.1, 9.3, 9.6 and 7, are **reproduced exactly**.

Unexpected discrepancies:
1. **PINN training depends on the CPU thread count.** With the baseline's
   2 threads the training is bit-reproducible. With 4 threads the
   floating-point reduction order changes and the optimization trajectory
   diverges. The resulting difference in final error (0.3 % relative) is
   far smaller than the network-seed spread in README §9.10 (0.487–0.567),
   so no PINN conclusion changes. The thread count is now recorded in
   every PINN run's metadata.
2. **Stored PINN evaluation metric, 7.5×10⁻⁷ difference.** The bit-identical
   model evaluates today to 0.48850211 for 1, 2 and 4 threads and on
   repeated calls. The stored baseline value is 0.48850136. The cause was
   not identified: the original evaluation ran under heavy CPU contention
   in an earlier session. The size is float32 round-off and is
   scientifically negligible.

Runtime and resources (4-core CPU, no GPU; `resources.jsonl`):

| Step | Wall time | CPU time (user) | Peak memory |
|---|---|---|---|
| validation | 11.4 min | 22.8 min | 0.67 GB |
| strength sweep (384 realizations) | 15.2 min | 31.0 min | 0.69 GB |
| correlation-length sweep (336) | 14.0 min | 28.9 min | 0.69 GB |
| localization (98 runs) | 1.9 min | 3.2 min | 0.70 GB |
| PINN, 4 threads | 4.7 min | 16.4 min | 1.95 GB |
| PINN, 2 threads | 7.8 min | 13.2 min | 1.97 GB |

Reference transport cost: ≈ 2.4 s wall per realization (L_d = 1, t_end = 1).

## 6. Remaining reproducibility risks

1. **Committed baseline results have no provenance metadata.** They are now
   treated as unverifiable caches, and the drivers refuse to reuse them. Four
   of them (validation, strength, correlation length, localization) and one
   PINN run now have bit-identical, fingerprinted twins in
   `results/reproduction/`. **Not yet reproduced:** the length, momentum,
   barrier and 192-realization convergence sweeps, and the other 17 PINN
   runs.
2. **Bit-level reproducibility was verified on one machine and one set of
   library versions** (`requirements-lock.txt`). Multithreaded FFTs and BLAS
   may differ in the last bits on other hardware or versions. For the
   reference solver this is far below the statistical errors (~10⁻²).
   PINN training can diverge more (item 3).
3. **PINN training is sensitive to the thread count** (Section 5). Exact PINN
   reproduction needs the same thread count. The count is recorded but is
   not part of the fingerprint.
4. **The disorder seed derivation** uses a constant salt in `disorder.py`
   (`SeedSequence([seed, 20241105])`). Changing it changes every
   realization. It is covered by the code fingerprint.
5. **Code changes only warn by default.** A cached result produced by
   different source code is reused with a warning unless `--strict-code` is
   given. Use `--strict-code` for any result that enters a report.
6. **Parameter lists of the sweeps live in the driver code.** They are
   covered by the code fingerprint and listed in each sweep's
   `run_metadata.json`, but they are not a separate configuration file.
7. `figures/pinn_vs_reference_W56.png` needs `eval_fields.npz`, which is
   git-ignored (it is regenerated by `evaluate.py` / `experiment_pinn.py`).
8. The reproduction runs record `git_dirty: true`: they ran before the
   provenance code was committed. The code fingerprint pins the exact
   sources, and the physics files are identical to `origin/main`.

## 7. Recommended next steps

1. Reproduce the remaining committed sweeps (length, momentum, barrier,
   convergence) into fresh directories with `--strict-code` (≈ 1.5–2 h of
   CPU). Then adopt the fingerprinted outputs as the canonical results.
2. Rerun the PINN study with a fixed, recorded thread count (e.g.
   `--threads 2`), and decide whether the thread count should become part of
   the PINN fingerprint.
3. Move the sweep parameter lists into versioned configuration files (e.g.
   JSON per sweep) so a sweep is fully described by data, not code.
4. Add a CI job that runs `python -m pytest -q tests` and a short
   reproduction (validation V1–V3 + one tiny transport point) against the
   lock file.
5. Only then proceed to Phase 2 (new physics: longer slabs, channel-number
   dependence, fixed-energy cross-checks), as listed in README §9.11.
