#!/usr/bin/env bash
# Reproduce the computer-assisted proof of the Auluck-Chowla-Gupta conjecture.
#   ./run_all.sh            analytic part + finite exact check to N=100000
#   ./run_all.sh --full     analytic part + finite exact check to N=200000 (the range used in the proof)
#   ./run_all.sh --quick    analytic part only (no finite check)
# Exits nonzero on any failure or mismatch against analytic/expected/ (timing tokens ignored).
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"
N=100000; QUICK=0
for a in "$@"; do case "$a" in
  --quick) QUICK=1;; --full) N=200000;; *) echo "unknown option $a"; exit 64;; esac; done
PY="${PYTHON:-python3}"
unset PYTHONOPTIMIZE   # the certificates rely on live assert statements
FAIL=0
norm() { # strip timing information: "time 12.3s" lines and "(12s)" / "(12.3s)" tokens
  sed -E -e '/^time [0-9.]+s$/d' -e 's/\(([0-9]+(\.[0-9]+)?)s\)/(Ts)/g' "$1"; }

# ---------- build ----------
make -s -C "$ROOT/finite" || { echo "BUILD FAILED"; exit 1; }

# ---------- analytic part (n >= 1e5) ----------
ACT="$ROOT/analytic/actual"; EXP="$ROOT/analytic/expected"; mkdir -p "$ACT"
cd "$ROOT/analytic"
run() { # name outfile args...
  local out="$1"; shift; local t0=$(date +%s)
  if ! "$PY" "$@" > "$ACT/$out" 2> "$ACT/$out.err"; then
    echo "FAIL (exit code)  $*"; FAIL=1; return; fi
  local dt=$(( $(date +%s) - t0 ))
  if diff <(norm "$EXP/$out") <(norm "$ACT/$out") > "$ACT/$out.diff"; then
    printf 'OK   %4ds  %s\n' "$dt" "$*"; rm -f "$ACT/$out.diff"
  else echo "MISMATCH  $*  (see analytic/actual/$out.diff)"; FAIL=1; fi
}
run ws_mech_elem.out    ws_mech_elem.py
run ws_mech_tv.out      ws_mech_tv.py
run ws_mech_lemA.out    ws_mech_lemA.py
run ws_mech_window.out  ws_mech_window.py 100000 5
run ws_mech_minor.out   ws_mech_minor.py 0.0042 3.0 0.4
run ws_p12_window.out   ws_p12_window.py 100000 5
run ws_p12_minor.out    ws_p12_minor.py 0.0042 3.0
run ws_p12_regionI.out  ws_p12_boxes.py 0.002 0.0040939 -2.52 2.10 3.1158 0.15926 1.01 0.01
run ws_p12_regionII.out ws_p12_asym.py 0.002 -2.52 2.10 0.15926 0.05
run ws_b_step.out       ws_b_step.py 0.002 0.0040939 -2.52 2.10 3.1158 1.02 0.02
run ws_b_asym.out       ws_b_asym.py 0.002 -2.52 2.10 0.05
run ws_b_edges.out      ws_b_edges.py 100000 0.002 0.0040939 0.14236
run ws_mech_edges.out   ws_mech_edges.py
run ws_cd_window.out    ws_cd_window.py 100000 5
run ws_cd_CR.out        ws_cd_CR.py 100000
run ws_cd_CLa.out       ws_cd_CLa.py 47 0.1
run ws_cd_CLb.out       ws_cd_CLb.py 100000 0.1 78
run ws_cd_D.out         ws_cd_D.py
run ws_mech_cd.out      ws_mech_cd.py 100000
run ws_mech_cla.out     ws_mech_cla.py 47 0.1
run ws_mech_mono.out    ws_mech_mono.py
# closing scripts for the formerly hand-proved items (notes/PROOFS_FULL.md)
run ws_close_rep.out    ws_close_rep.py
run ws_close_D.out      ws_close_D.py
run ws_close_elem.out   ws_close_elem.py
# threshold/hypothesis audit: reads the certified values from the outputs just produced in actual/
run ws_close_thresh.out ws_close_thresh.py "$ACT" ../notes/FiniteCheckExact.md
# audits (non-essential cross-checks)
run ws_b_minorM.out     ws_b_minorM.py 0.0042 3.0
run ws_b_check.out      ws_b_check.py 100000 1210 1358 1500

# ---------- finite part (n <= N) ----------
if [ "$QUICK" = 0 ]; then
  cd "$ROOT/finite"; mkdir -p out
  t0=$(date +%s)
  ./ws_fc_exact 2000 out/ex2000.txt && cmp out/ex2000.txt logs/ref2k.txt \
    && echo "OK   finite exact N=2000 == Python big-int reference" || { echo "FAIL finite N=2000"; FAIL=1; }
  if ./ws_fc_exact "$N" out/exlog.txt; then
    ref="logs/log$((N/1000))k.txt"
    if cmp out/exlog.txt "$ref"; then
      printf 'OK   %4ds  finite exact check N=%d: PASS, log byte-identical to %s\n' $(( $(date +%s)-t0 )) "$N" "$ref"
      sha256sum out/exlog.txt
    else echo "MISMATCH finite log vs $ref"; FAIL=1; fi
  else echo "FAIL finite exact check (exit $?)"; FAIL=1; fi
fi

if [ "$FAIL" = 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; fi
exit $FAIL
