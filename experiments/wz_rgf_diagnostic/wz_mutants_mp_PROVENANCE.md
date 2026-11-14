# Provenance of `wz_mutants_mp.py` (§4.3 transcription)

## What was supplied
The implementer received exactly one input file, `wz_mutants_mp_MATERIALS.md`
(sha256 `91fc23c523b1ea21914b2ffd3109dfe59dd156868ed77368ee09ea6c9202ca25`). It contains:

* the variable definitions;
* the §4.2 correct-code mathematical listing, verbatim;
* the N = 1 rules;
* the §4.2 **change column** only (18 rows; the oracle column was withheld);
* the required interface `mp_endpoint_blocks(Delta, lam_minus_1, mutant=None, dps=60)`, its return keys and trace
  names, and the precision rules (a private mpmath context; mpmath and stdlib only).

It did not receive `wz_mutants.py`, `wz_rgf.py`, `coupled_blocks.py`, `wz_mutants_oracle.py`, the T6 reference
(`cs2_t6_reference/`), or any test.

## Who wrote it, and how
`wz_mutants_mp.py` (sha256 at creation `a7339d25bdc1c621896fac6c2c23b4bed4cf2649bf1043f67f42b4c812fdfc8a`) was written by
a separately started AI subagent, in a fresh context that did not contain this session's conversation.

* **Instructions:** read only the materials file, and write scratch files only under a private scratch directory.
* **What it reported:** it read only the materials file and its own scratch files (`check.py`, `check2.py`). It checked
  its correct code against a dense mpmath inverse and found the corner blocks agree to 2.7e-60.
* **Interpretations it stated:**
  * the trace is recorded in actual computation order for M4/M4L/M4R/M9/M16;
  * M14's sign flip propagates literally through the chains;
  * `MUTANTS` follows the table order;
  * M12 differs from the correct code only by rounding;
  * at N_y = 1 the order mutants are value-equivalent.

The orchestrating session did not edit the file afterwards. Its sha256 above was recorded before any test was written.

## Independence: honest disclosure
Full independence **cannot be established**, for these reasons:

1. **Same model, same organisation.** The subagent is the same model as the author of `wz_mutants.py` and of the oracle.
   Shared priors mean shared misreadings of the spec are possible and are not excluded.
2. **Instruction-only isolation.** The subagent shared the container filesystem. Its restriction to the materials file
   was enforced by instruction, not by sandboxing, so its self-report cannot be verified externally.
3. **Spec-level coupling.** The materials file was prepared by the orchestrating session, which had read `wz_mutants.py`.
   It is a verbatim extract of the frozen spec plus an interface description, but its selection and wording were not
   produced independently.

## Mechanical guards
Enforced by `tests/test_wz_mutants_mp.py::test_transcription_and_oracle_import_independence`:

* the transcription imports only `mpmath` (AST), with no `__import__`, `importlib`, `exec` or `eval`;
* the transcription imports nothing from the solvers, the oracle or the T6 reference;
* the oracle does not import the solvers or the transcription.

## Agreement evidence
These are software checks, not scientific validation.

* `tests/test_wz_mutants_mp.py` (51 tests):
  * the correct code equals the exact direct inverse blocks within 1e-50·‖Ĝ‖;
  * every on-path mutant's designated block (plus M14's secondary blocks) agrees with the §4.2 oracle on exact
    T6 lattices X0/X1/X2/XD, all N in the T6 set, s ∈ {1, 1/2}, within 1e-50·‖Ĝ‖ entrywise;
  * the first deviating intermediate for M1 is `G_21^(2)`, equal to Ĝ^(1)_11 P_2;
  * determinism and private-context use; invalid labels are rejected.
* `tests/test_wz_mutants_oracle.py` (28 tests): the oracle agrees exactly with the exact rational recursion (rmutants).
