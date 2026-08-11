"""
describe19.py — the 19-vertex graph in closed form, and its certificate.

Paper A's G_17 is the three-hexagon union of Parshall/Mixon/Alexeev (21 points)
with four points deleted.  We found that putting two of them back raises the
certified bound on p^(4)_2 from 0.277518 to 0.288699, so the graph has to be
described explicitly rather than as "G_17 plus two vertices".

This script names every point symbolically in terms of

    w1 = e^{i pi/3} = 1/2 + i sqrt3/2        (the hexagonal rotation)
    w3 = 5/6 + i sqrt11/6                    (the Moser angle, cos = 5/6)

prints which four the 17-vertex version drops and which two the 19-vertex
version restores, and writes the graph plus certificate to JSON in the format
check_certificate.py expects.
"""
import sys, json, itertools
from fractions import Fraction as F
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog
from bigbuild2 import DEN, fadd, fsub, fmul, cmul_scaled
from build_variants import W1, W3, ONE, conj

TU=(DEN*DEN,0,0,0); T2=(4*DEN*DEN,0,0,0)
def cadd(a,b): return (fadd(a[0],b[0]), fadd(a[1],b[1]))
def csub(a,b): return (fsub(a[0],b[0]), fsub(a[1],b[1]))
def cn2(p): return fadd(fmul(p[0],p[0]), fmul(p[1],p[1]))
def key(p): return p[0]+p[1]

# ---- the 21 points, each with a symbolic name ---------------------------
w1c = conj(W1)
hexs = [(((0,)*4,(0,)*4), "0"),
        (ONE, "1"),
        (W1, "w1"),
        (cadd(ONE,W1), "1+w1"),
        (w1c, "conj(w1)"),
        (cadd(ONE,w1c), "1+conj(w1)"),
        (cadd(ONE,ONE), "2")]
v = cadd(ONE, W1)
named = list(hexs)
named += [(cmul_scaled(W3,p), f"w3*({nm})") for p,nm in hexs]
named += [(cadd(v, cmul_scaled(W3, csub(p,v))), f"v + w3*(({nm}) - v),  v = 1+w1")
          for p,nm in hexs]
named += [(cadd(cadd(ONE,ONE),W1), "2+w1"), (cadd(cadd(ONE,ONE),ONE), "3")]
uniq = {}
for p,nm in named:
    uniq.setdefault(key(p), (p,nm))
pts = list(uniq.values())
print(f"three-hexagon union: {len(pts)} distinct points "
      f"{'OK' if len(pts)==21 else 'MISMATCH'}")
ks = [key(p) for p,_ in pts]

def upairs(pl):
    return [(i,j) for i in range(len(pl)) for j in range(i+1,len(pl))
            if cn2(csub(pl[i],pl[j]))==TU]

def structure(pl):
    n=len(pl); units=[]; classes=defaultdict(list)
    for i in range(n):
        for j in range(i+1,n):
            v_=cn2(csub(pl[i],pl[j]))
            (units if v_==TU else classes[v_]).append((i,j))
    return units,classes

def profiles(n,units,pid,k=4):
    adj=[0]*n
    for u,v_ in units: adj[u]|=1<<v_; adj[v_]|=1<<u
    col=[-1]*n; out=set()
    def rec(i,used):
        if i==n:
            m=0
            for (a,b),t in pid.items():
                if col[a]==col[b]: m|=1<<t
            out.add(m); return
        forb=0
        for u in range(i):
            if adj[i]>>u&1: forb|=1<<col[u]
        for c in range(min(used+1,k)):
            if forb>>c&1: continue
            col[i]=c; rec(i+1,max(used,c+1)); col[i]=-1
    rec(0,0); return sorted(out)

def lp(idxs):
    pl=[pts[i][0] for i in idxs]
    units,classes=structure(pl)
    if T2 not in classes: return None
    pid={}
    for c,prs in classes.items():
        for p in prs: pid[p]=len(pid)
    profs=profiles(len(pl),units,pid)
    npair=len(pid); nv=npair+1
    cv=np.zeros(nv); cv[-1]=-1.0
    A=[];b=[]
    for m in profs:
        row=np.zeros(nv); mm=m; j=0
        while mm:
            if mm&1: row[j]=1.0
            mm>>=1; j+=1
        row[-1]=1.0; A.append(row); b.append(0.0)
    for c,prs in classes.items():
        row=np.zeros(nv)
        for p in prs: row[pid[p]]=-1.0
        A.append(row); b.append(1.0 if c==T2 else 0.0)
    r=linprog(cv,A_ub=np.array(A),b_ub=np.array(b),
              bounds=[(-10,10)]*npair+[(0,2)],method="highs")
    if not r.success: return None
    return r.x[-1], r.x[:npair], pid, classes, profs, units

# ---- recover G_17 by its signature --------------------------------------
g17=None
for drop in itertools.combinations(range(21),4):
    idxs=[i for i in range(21) if i not in drop]
    pl=[pts[i][0] for i in idxs]
    if len(upairs(pl))!=38: continue
    cls={cn2(csub(pl[i],pl[j])) for i in range(17) for j in range(i+1,17)
         if cn2(csub(pl[i],pl[j]))!=TU}
    if len(cls)==25: g17=(idxs,drop); break
idxs17, dropped = g17
print("\nthe four points G_17 deletes:")
for d in dropped: print(f"   {pts[d][1]}")

best=None
for add in itertools.combinations(dropped,2):
    idxs=sorted(idxs17+list(add))
    r=lp(idxs)
    if r and (best is None or r[0]>best[0][0]):
        best=(r,idxs,add)
(res, idxs19, added)=best
print("\nthe two the 19-vertex graph restores:")
for a in added: print(f"   {pts[a][1]}")
val,yf,pid,classes,profs,units=res
print(f"\n19-vertex graph: {len(idxs19)} vertices, {len(units)} unit edges, "
      f"{len(classes)} distance classes, {len(profs)} profiles")
print(f"  LP bound {val:.9f}   (G_17: 0.278044, certificate 0.277518)")
sys.stdout.flush()

# ---- rational certificate ----------------------------------------------
bestcert=None
for DENOM in (7680, 30720, 122880, 491520):
    y=[F(int(round(v*DENOM)),DENOM) for v in yf]
    for c,prs in classes.items():
        need=F(-1) if c==T2 else F(0)
        s=sum(y[pid[p]] for p in prs)
        if s<need:
            d=(need-s)/len(prs)
            for p in prs: y[pid[p]]+=d
    if not all(sum(y[pid[p]] for p in prs)>=(F(-1) if c==T2 else F(0))
               for c,prs in classes.items()): continue
    t=None
    for m in profs:
        s=F(0); mm=m; j=0
        while mm:
            if mm&1: s+=y[j]
            mm>>=1; j+=1
        t=-s if t is None else min(t,-s)
    if bestcert is None or t>bestcert[0]: bestcert=(t,y,DENOM)
t,y,DENOM=bestcert
print(f"\nCERTIFIED  p^(4)_2 >= {t.numerator}/{t.denominator} = {float(t):.9f}")
print(f"  vs Theorem 4.5 (11/40 = 0.275)     : "
      f"{'IMPROVED' if t>F(11,40) else 'no'}")
print(f"  vs shipped certificate 0.277518238 : "
      f"{'IMPROVED' if float(t)>0.277518238 else 'no'}")

out={
 "t": f"{t.numerator}/{t.denominator}",
 "n": len(idxs19),
 "unit_edges": len(units),
 "classes": len(classes),
 "profiles": len(profs),
 "vertices": [{"name": pts[i][1],
               "re": list(pts[i][0][0]), "im": list(pts[i][0][1])}
              for i in idxs19],
 "DEN": DEN,
 "basis": ["1","sqrt3","sqrt11","sqrt33"],
 "pairs": [[str(list(classes.keys()).index(c)) if False else "",
            p[0], p[1], f"{y[pid[p]].numerator}/{y[pid[p]].denominator}"]
           for c,prs in classes.items() for p in prs],
}
json.dump(out, open("certificate19.json","w"), indent=1)
print("\nwritten certificate19.json")
