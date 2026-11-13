"""T3 - lead-root enclosure (§2.4.2): exact inequalities decide; the mpmath dps-300 check is supplementary only."""
from fractions import Fraction as F
import mpmath

from . import inputs as I
from .result import Results


def run(prod):
    R = Results("T3 lead-root enclosure")
    ell = 256
    supp = []
    ctx = mpmath.MPContext()
    ctx.dps = 300
    for d in I.T3_DELTAS + I.T3_NY4_DELTAS:
        g = 1 - d
        enc = prod.lead_enclosure(g)
        lo, hi, mu, e = enc["lo"], enc["hi"], enc["mu"], enc["e"]
        xi = (1 - g * g) if abs(g) < 1 else (g * g - 1)
        tag = f"T3 delta={d} (gamma={g})"
        R.add(f"{tag}: lo^2 <= xi <= hi^2 (exact)", lo * lo <= xi <= hi * hi, "CERT-VIOLATION")
        R.add(f"{tag}: 2^-ell <= hi - lo <= 2^(1-ell)", F(1, 2**ell) <= hi - lo <= F(2, 2**ell), "CERT-VIOLATION")
        R.add(f"{tag}: e > 0", e > 0, "CERT-VIOLATION")
        if abs(g) < 1:
            R.add(f"{tag}: open branch, |lambda| = 1 bracketed: gamma^2 + lo^2 <= 1 <= gamma^2 + hi^2",
                  enc["branch"] == "open" and g * g + lo * lo <= 1 <= g * g + hi * hi, "CERT-VIOLATION")
            R.add(f"{tag}: open midpoint root gamma + i mu", enc["lam_mid"] == (g, mu), "CERT-VIOLATION")
        else:
            R.add(f"{tag}: closed branch, lambda in (0, 1): 0 < gamma - hi and gamma - lo < 1",
                  enc["branch"] == "closed" and 0 < g - hi and g - lo < 1, "CERT-VIOLATION")
            R.add(f"{tag}: closed midpoint root gamma - mu", enc["lam_mid"] == (g - mu, F(0)), "CERT-VIOLATION")
        R.add(f"{tag}: midpoint differs from the root (mu^2 != xi), so lambda~ != lambda and A_e != 0", mu * mu != xi,
              "CERT-VIOLATION")
        # supplementary (never decides): mpmath dps-300 sqrt(xi) lies in [lo, hi]
        sx = ctx.sqrt(ctx.mpf(xi.numerator) / xi.denominator)
        supp.append(ctx.mpf(lo.numerator) / lo.denominator <= sx <= ctx.mpf(hi.numerator) / hi.denominator)
    R.add(f"T3 supplementary mpmath dps-300 sqrt(xi) in [lo, hi] for all {len(supp)} cases (report only: {all(supp)})", True)
    return R
