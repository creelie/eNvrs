#!/usr/bin/env bash
# Builds the verifier against CAPD and runs the five verified cases of the paper (Table 1).
# Usage: cap/run_all.sh [path to capd-config]
# Each case prints PASS only if every inequality of Proposition 3.1 is verified.
set -euo pipefail
cd "$(dirname "$0")"
CAPD_CONFIG=${1:-${CAPD_CONFIG:-capd-config}}
# capd-config prepends its own pkgconfig directory to PKG_CONFIG_PATH using $PATH_SEPARATOR
export PATH_SEPARATOR=${PATH_SEPARATOR:-:}
g++ -O2 -std=c++17 verify_step.cpp -I. $($CAPD_CONFIG --cflags) $($CAPD_CONFIG --libs) -o verify_step
mkdir -p logs
# name            KB KA  Lambda  alpha_hi   centre                          radii                         SEXTRA
cases=(
 "L6_5_K80_7     80 7   6 5     377 1000   cases/center_L6_5_K80_7.txt     cases/radii_L6_5_K80_7.txt     0.02"
 "L13_10_K60_6   60 6   13 10   189 500    cases/center_L13_10_K60_6.txt   cases/radii_L13_10_K60_6.txt   0.02"
 "L3_2_K50_8     50 8   3 2     73 200     cases/center_L3_2_K50_8.txt     cases/radii_L3_2_K50_8.txt     0.05"
 "L17_10_K40_8   40 8   17 10   87 250     cases/center_L17_10_K40_8.txt   cases/radii_L17_10_K40_8.txt   0.02"
 "L9_5_K40_8     40 8   9 5     339 1000   cases/center_L9_5_K40_8.txt     cases/radii_L9_5_K40_8.txt     0.02"
)
fail=0
for c in "${cases[@]}"; do
  read -r name KB KA Ln Ld an ad cen rad sx <<<"$c"
  t0=$(date +%s.%N)
  if SEXTRA=$sx ./verify_step "$KB" "$KA" "$Ln" "$Ld" 0 1 "$an" "$ad" 1e-12 "$cen" "$rad" "logs/$name" 1.0 1e-31 > "logs/$name.out"; then st=PASS; else st=FAIL; fail=1; fi
  t1=$(date +%s.%N)
  printf '%-14s Lambda=%s/%s  alpha<=%s/%s  %s  (%.1f s)\n' "$name" "$Ln" "$Ld" "$an" "$ad" "$st" "$(echo "$t1 - $t0" | bc)"
  echo "seconds $(echo "$t1 - $t0" | bc)" >> "logs/$name.log"
done
exit $fail
