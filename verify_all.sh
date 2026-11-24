#!/usr/bin/env bash
# Re-check every computational claim in the repository.
#   ./verify_all.sh          fast checks (~3 min, plus ~2.5 min for chi>=7)
#   ./verify_all.sh --full   also the slow Paper B searches (~9 min more)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
run() { echo; echo "=== $1"; shift; (cd "$1" && shift && "$@"); }

run "Paper A: W6 colouring hull"            "$ROOT/code/paperA" python3 verify_hull.py
run "Paper A: G17 / G19 construction"       "$ROOT/code/paperA" python3 build_g17.py
run "Paper A: LP certificate (p5(sqrt3))"   "$ROOT/code/paperA" python3 check_certificate.py
run "Paper A: LP certificate (G19)"         "$ROOT/code/paperA" python3 check_certificate19.py
run "Paper B: n=8 excluded (h <= 7)"        "$ROOT/code/paperB" python3 h8_impossible.py
run "Paper B: exact 7-point configuration"  "$ROOT/code/paperB" python3 h7_theorem.py
run "Two-distance: upper bound 12"          "$ROOT/code/two_distance" python3 verify_upper12.py
run "Two-distance: G30 (chi = 6, critical)" "$ROOT/code/two_distance" python3 verify_g30.py
run "Two-distance: chi({1,d}) >= 7"         "$ROOT/code/two_distance" bash prove_chi7.sh

if [[ "${1:-}" == "--full" ]]; then
  run "Paper A: describe G19 distance classes" "$ROOT/code/paperA" env PYTHONPATH=../exploratory python3 describe19.py
  run "Paper B: sphere search"               "$ROOT/code/paperB" python3 h_sphere.py
  run "Paper B: n=7 exhaustive search"       "$ROOT/code/paperB" python3 h7.py
  run "Paper B: numerical Gram check"        "$ROOT/code/paperB" python3 h7_gram.py
fi
echo; echo "ALL CHECKS PASSED"
