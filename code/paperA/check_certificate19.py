"""
check_certificate19.py — standalone verification of the improved bound

        p^(4)_2  >=  2838001/9830400  =  0.2886964...

improving Theorem 4.5 (11/40 = 0.275) and the certificate shipped with it
(0.277518238).

The graph is the plain THREE-HEXAGON UNION of Parshall/Mixon/Alexeev:

    hexs  = { 0, 1, w1, 1+w1, conj(w1), 1+conj(w1), 2 }
    G_19  = hexs  u  w3*hexs  u  { v + w3*(p - v) : p in hexs },   v = 1 + w1

    w1 = e^{i pi/3} = 1/2 + i sqrt3/2 ,   w3 = 5/6 + i sqrt11/6

i.e. G_17 *without* Parshall's two extra points (2+w1 and 3) and *without*
Mixon's four deletions.  It has 19 vertices and 45 unit edges -- the counts that
build_g17.py already asserts on its first line.  The deletions that make G_17 a
minimal 4-chromatic graph turn out to cost the LP bound.

This script is independent of any LP solver.  It rebuilds the graph from the
definition in exact arithmetic over Q(sqrt3, sqrt11), enumerates every
monochromatic profile of a proper 4-colouring, reads certificate19.json, and
checks the two dual conditions in Fraction arithmetic:

    (C1)  -sum_{p monochromatic in profile} y_p  >=  t     for every profile
    (C2)  sum_{p in class c} y_p >= 0  for c != c*,   >= -1  for c = c*

Usage:  python3 check_certificate19.py
Needs:  certificate19.json in the same directory.  No floating point is used
        for any decision.
"""
import json, sys, itertools
from fractions import Fraction as F
from collections import defaultdict

# ---------- exact arithmetic in Q(sqrt3, sqrt11) --------------------------
# basis (1, sqrt3, sqrt11, sqrt33); every coordinate is an integer over DEN
DEN = 2**11 * 3**6
MUL = {}
for (i, j), (k, c) in {
    (0,0):(0,1), (0,1):(1,1), (0,2):(2,1), (0,3):(3,1),
    (1,1):(0,3), (1,2):(3,1), (1,3):(2,3),
    (2,2):(0,11), (2,3):(1,11), (3,3):(0,33),
}.items():
    MUL[(i,j)] = (k,c); MUL[(j,i)] = (k,c)

def fadd(x,y): return tuple(a+b for a,b in zip(x,y))
def fsub(x,y): return tuple(a-b for a,b in zip(x,y))
def fmul(x,y):
    o=[0,0,0,0]
    for i in range(4):
        if x[i]:
            for j in range(4):
                if y[j]:
                    k,c=MUL[(i,j)]; o[k]+=x[i]*y[j]*c
    return tuple(o)
def cadd(a,b): return (fadd(a[0],b[0]), fadd(a[1],b[1]))
def csub(a,b): return (fsub(a[0],b[0]), fsub(a[1],b[1]))
def cmul(a,b):
    re=fsub(fmul(a[0],b[0]), fmul(a[1],b[1]))
    im=fadd(fmul(a[0],b[1]), fmul(a[1],b[0]))
    assert all(v % DEN == 0 for v in re+im), "scale lost"
    return (tuple(v//DEN for v in re), tuple(v//DEN for v in im))
def cn2(p): return fadd(fmul(p[0],p[0]), fmul(p[1],p[1]))
def key(p): return p[0]+p[1]

ZERO=(0,0,0,0)
ONE=((DEN,0,0,0), ZERO)
W1=((DEN//2,0,0,0),(0,DEN//2,0,0))          # 1/2 + i sqrt3/2
W3=((5*DEN//6,0,0,0),(0,0,DEN//6,0))        # 5/6 + i sqrt11/6
W1C=(W1[0], tuple(-v for v in W1[1]))
UNIT=(DEN*DEN,0,0,0)
D2_TARGET=(4*DEN*DEN,0,0,0)                 # squared distance 4, i.e. |d| = 2

assert cn2(W1)==UNIT and cn2(W3)==UNIT, "w1 and w3 must be units"

# ---------- 1. the graph ---------------------------------------------------
hexs=[(ZERO,ZERO), ONE, W1, cadd(ONE,W1), W1C, cadd(ONE,W1C), cadd(ONE,ONE)]
v=cadd(ONE,W1)
allp = hexs + [cmul(W3,p) for p in hexs] \
            + [cadd(v, cmul(W3, csub(p,v))) for p in hexs]
uniq={}
for p in allp: uniq.setdefault(key(p), p)
pts=list(uniq.values())
n=len(pts)
units=[]; classes=defaultdict(list)
for i in range(n):
    for j in range(i+1,n):
        d=cn2(csub(pts[i],pts[j]))
        (units if d==UNIT else classes[d]).append((i,j))
print(f"graph: n={n}, unit edges={len(units)}, non-unit pairs="
      f"{sum(len(x) for x in classes.values())}, distance classes={len(classes)}")
assert n==19 and len(units)==45, "must be the 19-vertex / 45-edge union"
assert D2_TARGET in classes, "the distance-2 class must be present"

# no two classes may coincide numerically
from math import sqrt
def approx(d): return (d[0]+d[1]*sqrt(3)+d[2]*sqrt(11)+d[3]*sqrt(33))/(DEN*DEN)
vals=sorted((approx(d),d) for d in classes)
assert all(abs(vals[i][0]-vals[i+1][0])>1e-9 for i in range(len(vals)-1)), \
    "two distance classes coincide numerically"

# ---------- 2. all monochromatic profiles of proper 4-colourings -----------
pid={}
for c,prs in classes.items():
    for p in prs: pid[p]=len(pid)
adj=[0]*n
for u,w in units: adj[u]|=1<<w; adj[w]|=1<<u
col=[-1]*n; profs=set()
def rec(i,used):
    if i==n:
        m=0
        for (a,b),t in pid.items():
            if col[a]==col[b]: m|=1<<t
        profs.add(m); return
    forb=0
    for u in range(i):
        if adj[i]>>u&1: forb|=1<<col[u]
    for c in range(min(used+1,4)):
        if forb>>c&1: continue
        col[i]=c; rec(i+1,max(used,c+1)); col[i]=-1
rec(0,0)
profs=sorted(profs)
print(f"distinct monochromatic profiles of proper 4-colourings: {len(profs)}")

# ---------- 3. read and check the certificate -----------------------------
cert=json.load(open('certificate19.json'))
t=F(cert['t'])
y=[F(0)]*len(pid)
for _, a, b, val in cert['pairs']:
    y[pid[(a,b)]]=F(val)

ok2=True
for c,prs in classes.items():
    s=sum(y[pid[p]] for p in prs)
    need=F(-1) if c==D2_TARGET else F(0)
    if s<need:
        ok2=False
        print(f"  (C2) FAILS at |d|^2 = {approx(c):.6f}: sum y = {s} < {need}")
print(f"(C2)  {len(classes)} class conditions: {'PASS' if ok2 else 'FAIL'}")

worst=None
for m in profs:
    s=F(0); mm=m; j=0
    while mm:
        if mm&1: s+=y[j]
        mm>>=1; j+=1
    worst=-s if worst is None else min(worst,-s)
ok1=(t<=worst)
print(f"(C1)  {len(profs)} profile conditions: {'PASS' if ok1 else 'FAIL'}   "
      f"(min over profiles of -sum y = {float(worst):.9f})")

if ok1 and ok2:
    print(f"\nCERTIFIED:  p^(4)_2 >= {t} = {float(t):.9f}")
    print(f"            Theorem 4.5  11/40 = 0.275         -> improvement: {t>F(11,40)}")
    print(f"            shipped certificate 0.277518238    -> improvement: "
          f"{float(t)>0.277518238}")
    print(f"            published bound 1/4 = 0.25         -> improvement: {t>F(1,4)}")
else:
    print("\nCERTIFICATE REJECTED")
    sys.exit(1)
