"""POST-HOC (written after amplitude_roundoff.py ran): the relative r error of the production N = 2 device is ill-posed because
the exact r of its represented lattice (R_code, 1-ulp lead/device mismatch) is itself at rounding level. Report ABSOLUTE
|r_prod - r_ref| and |r_ref| for the open channels. Read-only use of stationary / validate_rgf_b1."""
import json
import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
from stationary import build_device, scattering_matrix  # noqa: E402
from amplitude_roundoff import exact_hat, mpc  # noqa: E402

mp.mp.dps = 110
E, dx, ny = 125.69231487269683, 1 / 1024, 63
dev2 = build_device(lambda X, Y: np.zeros_like(X), 1.7, 1.7 + dx, dx, 1.0, ny)
r2 = scattering_matrix(E, dev2)
t = float(dev2.t_x)
Ep = dev2.basis.E_perp
out = {}
for kk, n in enumerate(np.where(r2["modes"]["open"])[0]):
    a_dev = E - (2.0 * t + Ep[n])
    Rc = exact_hat(mp.mpf(a_dev) / (2 * t), -mpc(r2["modes"]["sigma"][n]) / t, 2)
    gam = mp.mpf(float(r2["modes"]["gamma"][n]))
    r_ref = 1j * gam * Rc["G11"] / t - 1
    out[int(n) + 1] = dict(abs_r_ref=float(abs(r_ref)), abs_r_prod=float(abs(r2["r_L"][kk, kk])),
                           abs_err_r=float(abs(mpc(r2["r_L"][kk, kk]) - r_ref)))
print(json.dumps(out, indent=1))
with open(os.path.join(HERE, "results_amplitude_roundoff", "posthoc_N2_reflection.json"), "w") as fh:
    json.dump(out, fh, indent=1)
