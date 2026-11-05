#!/bin/sh
# Baseline reproduction (Phase 1). Re-runs a representative subset of the
# committed baseline experiments into fresh directories under
# results/reproduction/ and records wall time and peak memory per step.
set -e
R=results/reproduction
step() {  # name, command...
  name=$1; shift
  python - "$name" "$@" <<'PY'
import resource, subprocess, sys, time, json
name, cmd = sys.argv[1], sys.argv[2:]
t0 = time.time()
p = subprocess.run(cmd)
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
rec = dict(step=name, cmd=" ".join(cmd), returncode=p.returncode, wall_s=round(time.time() - t0, 1),
           cpu_user_s=round(ru.ru_utime, 1), cpu_sys_s=round(ru.ru_stime, 1), peak_rss_mb=round(ru.ru_maxrss / 1024, 1))
open("results/reproduction/resources.jsonl", "a").write(json.dumps(rec) + "\n")
print("RESOURCE", json.dumps(rec), flush=True)
sys.exit(p.returncode)
PY
}
step validation   python validate_solver.py --out $R/validation
step strength     python experiment_disorder.py --sweep strength --n-real 48 --out $R/transport
step corr_length  python experiment_disorder.py --sweep corr_length --n-real 48 --out $R/transport
step localization python experiment_localization.py --n-real 16 --out $R/localization
step pinn_clean   python experiment_pinn.py --only clean_fourier_uniform_s0 --out $R/pinn
echo REPRODUCTION_DONE
