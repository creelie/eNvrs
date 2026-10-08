#!/usr/bin/env bash
# fast_checks.sh -- the independent checks that do not need CAPD (a few seconds each).
#
#   1. crosscheck/recheck.py   re-checks, in 60-digit interval arithmetic, the scalar
#                              inequalities of Proposition 3.1 from the logged enclosures
#                              in cap/logs/ and the input files in cap/cases/;
#   2. crosscheck/front_steps.c iterates the full model (all deep shells, 14 ahead shells,
#                              viscosity) for 150 steps from two corners of the box for
#                              Lambda = 13/10, alpha = 3/8 and requires the state to stay
#                              inside the box (non-rigorous, long double);
#   3. crosscheck/kappa.jl     recomputes kappa(Lambda) in Julia, if julia is installed,
#                              and compares it with numerics/kappa_table.txt.
# Needs python3 with mpmath, a C compiler; optionally julia.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/build/fast-checks"
mkdir -p "$OUT"
cd "$ROOT/crosscheck"

echo "== interval re-check of the scalar inequalities"
python3 recheck.py | tee "$OUT/recheck.log"
grep -q "^all scalar inequalities re-checked$" "$OUT/recheck.log"

echo "== full-model iteration in long double (Lambda = 13/10, alpha = 3/8, 150 steps)"
cc -O2 -o "$OUT/front_steps" front_steps.c -lm
for seed in 0 1; do
  "$OUT/front_steps" 13 10 0.375 1e-12 150 60 6 \
    ../cap/cases/center_L13_10_K60_6.txt ../cap/cases/radii_L13_10_K60_6.txt "$seed" \
    | tail -1 | tee -a "$OUT/front_steps.log"
done
awk '/worst ratio/ { for (i = 1; i <= NF; ++i) if ($i == "ratio") w = $(i + 1);
                     if (!(w < 1)) bad = 1 }
     END { if (bad) { print "front left the box" > "/dev/stderr"; exit 1 } }' "$OUT/front_steps.log"

if command -v julia > /dev/null; then
  echo "== kappa(Lambda) in Julia"
  julia kappa.jl | tee "$OUT/kappa.log"
  python3 - "$OUT/kappa.log" ../numerics/kappa_table.txt <<'PY'
import re, sys
tab = {}
for line in open(sys.argv[2]):
    if line.strip() and not line.startswith('#'):
        f = line.split(); tab[round(float(f[0]), 4)] = float(f[1])
worst = 0.0
for line in open(sys.argv[1]):
    m = re.match(r"Lambda=([0-9.]+) kappa=([0-9.]+)", line)
    if m:
        worst = max(worst, abs(float(m.group(2)) - tab[round(float(m.group(1)), 4)]))
print(f"largest difference from numerics/kappa_table.txt: {worst:.2e}")
sys.exit(0 if worst < 1e-9 else 1)
PY
else
  echo "== julia not installed; skipping the kappa check"
fi
echo "fast checks passed"
