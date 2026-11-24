#!/usr/bin/env bash
# Refute conf7.cnf with CaDiCaL (via python-sat), write a DRUP proof, and check it.
# Uses drat-trim if it is on PATH, and always the bundled independent checker rup_checker.c.
set -euo pipefail
cd "$(dirname "$0")"
python3 verify_chi7.py
python3 - <<'EOF'
from pysat.formula import CNF
from pysat.solvers import Solver
F = CNF(from_file="conf7.cnf")
s = Solver(name="cadical195", bootstrap_with=F.clauses, with_proof=True)
assert s.solve() is False, "expected UNSAT"
open("conf7.drup", "w").write("\n".join(s.get_proof()) + "\n")
print("CaDiCaL: UNSAT, proof written to conf7.drup")
EOF
gcc -O2 -o rup_checker rup_checker.c
out=$(./rup_checker conf7.cnf conf7.drup); echo "$out" | tail -2
echo "$out" | grep -q "VERIFIED UNSAT" || { echo "PROOF CHECK FAILED"; exit 1; }
if command -v drat-trim >/dev/null 2>&1; then drat-trim conf7.cnf conf7.drup | grep -E "^s "; fi
