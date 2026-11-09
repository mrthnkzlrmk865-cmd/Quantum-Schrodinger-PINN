"""
POST-HOC check (written AFTER run_worst_point.py was executed; not part of PREDICTIONS.md).

Question raised by the run: why is the W/Z (delta-lead) amplitude error eps_A = +4.82e-14 about 44x larger than
the original solver's (-1.08e-15) at the worst B1 point?
Hypothesis checked here: with a clean channel, W_j stays at one double value, so the per-slice propagator
P = fl((1 - W)^{-1}) is the same double at every slice; its modulus differs from 1 by the rounding of Re P ~ 1,
and Ghat_N1 ~ P^(N-1) turns that into a coherent amplitude drift (N-1)(|P|-1).
Arithmetic only on the double W/Z inputs; numpy + mpmath; no solver import.
"""
import json
import pathlib
import sys

import mpmath as mp
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wz_rgf import lead_from_delta  # noqa: E402

mp.mp.dps = 60
N = 3483
d = 2.21467958362565e-06                      # delta_fl (exact; see results/worst_point.json)
lm1 = lead_from_delta(np.array([d]))["lam_minus_1"][0]
W = 2.0 * d + lm1                             # W_1 in double
I1 = np.eye(1, dtype=complex)
P = np.linalg.solve(I1 - np.array([[W]]), I1)[0, 0]   # same operation as wz_endpoint_blocks

D = mp.mpf(d)
w_exact = mp.mpc(D, mp.sqrt(D * (2 - D)))     # exact w* = 1 - 1/lambda = delta + i sin(theta)
Wm = mp.mpc(W.real, W.imag)
Pm = mp.mpc(P.real, P.imag)
out = {
    "W_minus_w_exact_abs": float(abs(Wm - w_exact)),
    "abs_exact_inverse_of_double_W_minus_1": float(abs(1 / (1 - Wm)) - 1),
    "abs_P_double_minus_1": float(abs(Pm) - 1),
    "predicted_eps_A_coherent": float((N - 1) * (abs(Pm) - 1)),
    "observed_eps_A_WZ_delta_lead": json.load(open(pathlib.Path(__file__).parent / "results/worst_point.json"))
    ["comparisons"][2]["eps_A_tRL"],
    "arg_P_error": float(mp.arg(Pm) - mp.arg(1 / (1 - w_exact))),
}
out["predicted_eps_phi_coherent"] = (N - 1) * out["arg_P_error"]
print(json.dumps(out, indent=1))
(pathlib.Path(__file__).parent / "results/posthoc_amplitude.json").write_text(json.dumps(out, indent=1) + "\n")
