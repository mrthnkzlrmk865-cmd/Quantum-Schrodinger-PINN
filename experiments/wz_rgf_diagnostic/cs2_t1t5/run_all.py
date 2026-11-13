"""Stage 4: T1-T5 in order (§2.4.2)."""
from . import t1, t2, t3, t4, t5


def run(prod, cov):
    return [t1.run(prod), t2.run(prod, cov), t3.run(prod), t4.run(prod, cov), t5.run(prod)]
