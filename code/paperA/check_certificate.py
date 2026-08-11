"""
check_certificate.py — standalone verification of Theorem 4.5:  p^(4)_2 >= 11/40.

Rebuilds G_17 from its definition in Q(sqrt3, sqrt11), enumerates all monochromatic
profiles of proper 4-colourings, reads certificate.json, and checks the two dual
conditions in exact rational arithmetic (Fraction only, no floating point).

Usage:  python3 check_certificate.py
Needs:  sympy, and certificate.json in the same directory.
"""
import json, itertools
from fractions import Fraction as F
import sympy as sp

# ---------- 1. the graph -------------------------------------------------
I = sp.I
w1 = (1 + I*sp.sqrt(3))/2
w3 = (sp.Integer(5) + I*sp.sqrt(11))/6
hexs = [sp.Integer(0), sp.Integer(1), w1, 1+w1, sp.conjugate(w1), 1+sp.conjugate(w1), sp.Integer(2)]
v = 1 + w1
pts21 = hexs + [w3*p for p in hexs] + [v + w3*(p-v) for p in hexs] + [2+w1, sp.Integer(3)]

uniq = []
for p in pts21:
    p = sp.expand(sp.simplify(p))
    if not any(sp.simplify(p-q) == 0 for q in uniq):
        uniq.append(p)
assert len(uniq) == 21, len(uniq)

rem = [sp.Integer(3),
       sp.Rational(1, 2)*I*sp.sqrt(3) + sp.Rational(5, 2),
       sp.Rational(1, 12)*(I*sp.sqrt(11) - 1)*(-I*sp.sqrt(3) - 3),
       sp.Rational(1, 12)*(-I*sp.sqrt(11) - 5)*(-I*sp.sqrt(3) - 1)
       + sp.Rational(1, 12)*(I*sp.sqrt(11) - 1)*(-I*sp.sqrt(3) - 3)]
rem = [sp.expand(sp.simplify(r)) for r in rem]
pts = [p for p in uniq if not any(sp.simplify(p-r) == 0 for r in rem)]
n = len(pts)
assert n == 17, n

d2 = {}
for a, b in itertools.combinations(range(n), 2):
    d2[(a, b)] = sp.simplify(sp.expand(sp.Abs(pts[a]-pts[b])**2))
edges = [e for e, d in d2.items() if d == 1]
assert len(edges) == 38, len(edges)

classes = {}
for e, d in d2.items():
    if d != 1:
        classes.setdefault(sp.srepr(d), []).append(e)
names = sorted(classes)
cstar = [k for k in names if sp.sympify(k) == 4][0]
print(f"graph: n={n}, unit edges={len(edges)}, non-unit pairs={sum(map(len,classes.values()))}, "
      f"distance classes={len(names)}")

# no two classes may coincide numerically
vals = sorted((float(sp.N(sp.sympify(k))), k) for k in names)
assert all(abs(vals[i][0]-vals[i+1][0]) > 1e-12 for i in range(len(vals)-1)), "class collision"

# ---------- 2. all monochromatic profiles of proper 4-colourings ---------
adj = [0]*n
for a, b in edges:
    adj[a] |= 1 << b; adj[b] |= 1 << a
order = sorted(range(n), key=lambda x: -bin(adj[x]).count('1'))
pos = {x: i for i, x in enumerate(order)}
prev = [[u for u in range(n) if (adj[x] >> u) & 1 and pos[u] < pos[x]] for x in range(n)]
pairs = [e for c in names for e in classes[c]]
idx = {e: i for i, e in enumerate(pairs)}
col = [-1]*n
profs = set()

def rec(i, used):
    if i == n:
        p = 0
        for e, j in idx.items():
            if col[e[0]] == col[e[1]]: p |= 1 << j
        profs.add(p); return
    x = order[i]; forb = 0
    for u in prev[x]: forb |= 1 << col[u]
    for c in range(min(used+1, 4)):
        if (forb >> c) & 1: continue
        col[x] = c; rec(i+1, max(used, c+1))
    col[x] = -1

rec(0, 0)
profs = sorted(profs)
print(f"distinct monochromatic profiles of proper 4-colourings: {len(profs)}")

# ---------- 3. read and check the certificate ---------------------------
cert = json.load(open('certificate.json'))
t = F(cert['t'])
y = [F(0)]*len(pairs)
for dstr, a, b, val in cert['pairs']:
    y[idx[(a, b)]] = F(val)

ok2 = True
for nm in names:
    s = sum(y[idx[e]] for e in classes[nm])
    need = F(-1) if nm == cstar else F(0)
    if s < need:
        ok2 = False
        print(f"  (C2) FAILS at d^2 = {sp.sympify(nm)}: sum y = {s} < {need}")
print(f"(C2)  25 class conditions: {'PASS' if ok2 else 'FAIL'}")

worst = None
for p in profs:
    s = F(0); pp = p; j = 0
    while pp:
        if pp & 1: s += y[j]
        pp >>= 1; j += 1
    worst = -s if worst is None else min(worst, -s)
ok1 = (t <= worst)
print(f"(C1)  {len(profs)} profile conditions: {'PASS' if ok1 else 'FAIL'}   "
      f"(min over profiles of -sum y = {float(worst):.9f})")

if ok1 and ok2:
    print(f"\nCERTIFIED:  p^(4)_2 >= {t} = {float(t):.9f}")
    print(f"            11/40 = 0.275   -> Theorem 4.5 {'holds' if t >= F(11,40) else 'FAILS'}")
    print(f"            published bound 1/4 = 0.25  -> improvement: {t > F(1,4)}")
