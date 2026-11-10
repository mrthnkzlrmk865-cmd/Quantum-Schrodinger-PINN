# Predictions: W/Z amplitude-roundoff investigation (written and committed BEFORE any of its code was run)

Production baseline: `74227189c524fb346247c92dc0a9602acf64730d`. Prior experiment: `edf3881` (worst point, scalar 1x1 path only;
see `ERRATUM_worst_point.md`). Point unchanged: B1 config F, N_y = 63, dx = 1/1024, N = 3483, E = 125.69231487269683.
Inputs inspected (not a diagnostic run): the F device has V_j = 0 exactly (max |V| = 0.0), 5 open channels (n = 1..5),
58 closed; delta_n = (E - E_perp,n) dx^2 ranges from +1.15e-4 (n = 1) through +2.21e-6 (n = 5) to -4.96e-5 (n = 6, closed,
nearest threshold) and -1.856e-2 (n = 63, deepest closed). The numbers below use only arithmetic on these inputs, the
committed results of `edf3881`, and the post-hoc check already published there (which is therefore NOT independent evidence).

## Hypothesis under test (H_A)
The W/Z amplitude error at the worst point is a coherent drift caused by the modulus defect of the rounded per-slice
propagator: with a clean channel, W_j is one fixed double, P_d = fl((1 - W)^-1) is one fixed double with |P_d| != |lambda_W| = 1,
and G_hat_{N-1,1} contains N - 1 such factors, so eps_A ~= sum_j (|P_d,j| / |lambda_W| - 1) ~= (N-1) d_A.

Pre-registered test (exact complex-log telescoping of t_RL, residual computed at >= 100 digits):
ln(t_d / t_ref) = [W recurrence] + [P formation] + [propagation multiplications] + [Z_N formation] + [Z_N inversion]
+ [lead data in Fisher-Lee] + [Fisher-Lee multiplication]. Real parts = amplitude contributions, imaginary parts = phase.
- H_A SUPPORTED if, for every inversion method on channel 5 (1x1 and 63x63 paths), the "P formation" real term explains the
  final ln|t_d/t_ref| to within 10 % + 2e-15, and the W-recurrence term is < 10 % of it.
- H_A REFUTED if, for any method, another stage contributes >= 50 % of the amplitude error.
- Otherwise PARTIAL. The test is informative only if at least two methods produce distinct P_d; if all methods give
  bit-identical P_d the result will be reported as "consistent but not independently discriminated".

## Expected source of amplitude drift (channel 5)
- P formation: d_A = |P_d| - 1 of order 1e-17 (|1/(1 - W_d)| - 1 is of order 1e-22, i.e. negligible). Re P ~ 1 is
  representable only to ulp(1)/2 ~ 1.1e-16, so for ANY double-valued P the modulus defect is bounded by ~1.1e-16 per step,
  and the coherent worst case is (N-1) * 1.1e-16 ~ 3.9e-13. Expected P-formation term for np.linalg.solve: +4.76e-14
  (= the post-hoc value; reproduction, not independent).
- W recurrence: |W_j - w*| ~ 1e-19 and constant in j; contribution to ln|t| <~ 1e-15.
- Propagation multiplications (P @ G, N - 2 of them): the multiplicand rotates, so rounding is not coherent;
  expected |contribution| ~ sqrt(N) * 1e-16 ~ 6e-15 or less.
- Z_N formation and inversion, Fisher-Lee lead data and multiplication: O(1e-16..1e-15) each (single operations;
  Z_N has |Z_N| = 2 sin theta = 4.2e-3 but relative rounding stays O(eps)).

## Expected difference between 1x1 and 63x63 solves
The 63x63 clean Delta_j is exactly diagonal, so LU needs no pivoting and every BLAS/LAPACK product involves only exact zeros
off the diagonal. Prediction: off-diagonal entries of W_j, P_j, G_hat and Z_N^-1 stay exactly 0. Channel-5 results are
expected (most likely) to be bit-identical between the 1x1 and 63x63 paths. This is uncertain, because LAPACK/OpenBLAS may
use size-dependent kernels (e.g. reciprocal-then-multiply in trsm, FMA in gemm). If they differ, the difference is
expected at per-step rounding level, giving |eps_A(63) - eps_A(1)| <= (N-1) * 1.1e-16 ~ 3.9e-13.
A diagonal clean matrix cannot validate the multiplication order needed for coupled (disordered) blocks; that remains untested.

## Expected effect of the inversion method and of iterative refinement
- solve, LU (scipy lu_factor/lu_solve), inv: expected to give P_d within 1 ulp per component of each other; they may or may
  not be bit-identical (no prior assumption about which is better).
- Scalar complex division (diagonal/1x1 control only): expected within 1 ulp of solve.
- One step of iterative refinement with the residual R = I - A P computed in double: expected to change P_d by at most
  ~1 ulp per component and NOT to reduce |eps_A| systematically (the residual is itself at rounding level).
- Refinement with the residual computed exactly (mpmath) and rounded to double (1x1 only): expected to give (close to)
  the correctly rounded P. Because the correctly rounded P also has |P| != 1 at the ~1e-17..1e-16 level, it is expected
  NOT to remove the drift. It changes eps_A only if solve's P is not already correctly rounded.
- No method is expected to bring |eps_A| for channel 5 below ~1e-15 unless its P_d happens to have |P_d| - 1 below ~3e-19.
- No normalization of |P| is applied anywhere.

## Expected effect on phase accuracy
Per-step phase defect of P_d is <~ ulp(Im P)/2 + Im P * ulp(1)/2 ~ 5e-19. Final |eps_phi| for channel 5 is expected
<= 2e-15 for every method and both paths (prior 1x1 solve: -1.19e-15). The open channels n = 1..4 are expected at
|eps_phi| <= 2e-15 as well. The 1/sin(theta) amplification of the original formulation is absent here.

## Expected unitarity and reciprocity
- Full open-channel S (10x10, 5 open channels), entrywise max of |S^dag S - I| and |S S^dag - I|: dominated by
  2 max_n |eps_A,n|. Expected ~1e-13 (<= 8e-13) for W/Z; the production solver's full S is expected ~1e-13 as well.
- Reciprocity max |S - S^T|: expected exactly 0 for W/Z in the clean diagonal case (left and right chains multiply the same
  diagonal numbers; complex multiplication is commutative in IEEE arithmetic). Production: expected ~1e-14 (prior 2x2: 5.9e-15).
- Off-diagonal entries of the W/Z open-channel S: exactly 0. Off-diagonal entries of the production S: at rounding level
  or exactly 0 (prior max = 0.0).

## Closed channels and conditioning (63x63 path)
- Closed lead root, delta-based formula: relative error of (lambda - 1) <= ~4 eps against the exact root of the
  represented delta. The production root c - sqrt(c^2 - 1) suffers cancellation in c^2 - 1: relative error of (lambda - 1)
  up to ~eps/|delta| (~2e-12 for n = 6, ~1e-14 for n = 63) against the exact root of its own represented c.
- Closed-channel G_hat blocks: relative error of G_hat_11 and G_hat_NN <= 1e-14. Relative error of G_hat_N1 (a product of
  N - 1 contractions) <= N eps ~ 8e-13.
- max_j |W_j| (entrywise) = |w*_63| = |delta_63| + s_63 ~= 0.2121.
- cond_2(I - W_j) ~= max_n |1/lambda_n| / min_n |1/lambda_n| = 1.2121 / 1 ~= 1.212 at every slice.
- cond_2(Z_N) ~= max_n |lambda_n - 1/lambda_n| / min_n |...| = 0.3871 / 0.0042092 ~= 92.

## Edge cases (Part D): N = 1 and N = 2, all 63 channels
- High-precision W/Z equations vs high-precision original recursion vs direct inversion of M vs closed form on the same
  lattice: agreement to <= 1e-90 (algebraic identity).
- Double W/Z vs exact (its own represented lattice): relative error <= 1e-14 for every channel (open, closed near
  threshold n = 6, closed far n = 63). Open S: U <= 1e-14, reciprocity exactly 0.
- Production original at N = 2 vs its own represented lattice: relative error <= 1e-13 for the open channels. Production cannot
  build an N = 1 device (build_device requires >= 1 interval), so N = 1 is tested for W/Z and the mp formulations only.

## Interpretation rule
Whatever the outcome: one energy, one device, clean only; no validation claim; B1 still FAILS at F under the unchanged 1e-10.
