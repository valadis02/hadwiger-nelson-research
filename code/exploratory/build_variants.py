"""
build_variants.py — three containers of comparable size, three DIFFERENT
generating strategies, so that the exploration has something to compare against.

  A  TRANSLATION closure, Moser angle.       dirs = w1^a w3^b, |b|<=2   (ours)
  B  TRANSLATION closure, de Grey's angle.   dirs = w1^a w4^b, |b|<=2
       w4 = e^{i*arcsin(1/sqrt12)} : cos = sqrt33/6, sin = sqrt3/6.
       Same field Q(sqrt3,sqrt11) since sqrt33 = sqrt3*sqrt11, so the same
       integer machinery applies unchanged.
  C  ROTATION closure (spindling).           orbit of a seed under rotations
       about its OWN vertices.  Structurally different: A and B are closed
       under translation by units, C is not; C is closed under rotation, which
       is the operation that built the Moser spindle in the first place.

Everything exact, integers only.  DEN must clear every angle used.
"""
import time, sys, resource
from bigbuild2 import fadd, fsub, fmul, cmul_scaled, DEN

def mkpt(re, im): return (re, im)

ONE = ((DEN,0,0,0),(0,0,0,0))
W1  = ((DEN//2,0,0,0),(0,DEN//2,0,0))                 # 1/2 + i sqrt3/2
W3  = ((5*DEN//6,0,0,0),(0,0,DEN//6,0))               # 5/6 + i sqrt11/6   (Moser)
W4  = ((0,0,0,DEN//6),(0,DEN//6,0,0))                 # sqrt33/6 + i sqrt3/6 (de Grey)

def conj(z): return (z[0], tuple(-v for v in z[1]))

def check_unit(z, name):
    n2 = fadd(fmul(z[0],z[0]), fmul(z[1],z[1]))
    assert n2 == (DEN*DEN,0,0,0), (name, n2)
    print(f"   |{name}| = 1 exactly  OK")

check_unit(W1,"w1"); check_unit(W3,"w3"); check_unit(W4,"w4")

def dirs_from(w, bmax):
    out={}
    wc=conj(w)
    for a in range(6):
        wa=ONE
        for _ in range(a): wa=cmul_scaled(wa,W1)
        for b in range(-bmax,bmax+1):
            z=wa
            for _ in range(abs(b)): z=cmul_scaled(z, w if b>0 else wc)
            out.setdefault(z[0]+z[1], z)
    return list(out.values())

def translation_closure(dirs, depth, cap):
    cur={(0,)*8: ((0,0,0,0),(0,0,0,0))}
    for s in range(depth):
        nxt={}
        for (re,im) in cur.values():
            for (dre,dim) in dirs:
                r=(re[0]+dre[0],re[1]+dre[1],re[2]+dre[2],re[3]+dre[3])
                m=(im[0]+dim[0],im[1]+dim[1],im[2]+dim[2],im[3]+dim[3])
                nxt[r+m]=(r,m)
        cur=nxt
        if len(cur)>cap: break
    return cur

def rotate_about(p, c, w):
    """c + w*(p-c), all scaled by DEN"""
    d=(fsub(p[0],c[0]), fsub(p[1],c[1]))
    rd=cmul_scaled(d,w)
    return (fadd(c[0],rd[0]), fadd(c[1],rd[1]))

def rotation_closure(seed, angles, rounds, cap):
    cur={p[0]+p[1]: p for p in seed}
    for r in range(rounds):
        centres=list(cur.values())
        nxt=dict(cur)
        for c in centres:
            for w in angles:
                for p in centres:
                    q=rotate_about(p,c,w)
                    nxt[q[0]+q[1]]=q
                    if len(nxt)>cap: break
                if len(nxt)>cap: break
            if len(nxt)>cap: break
        cur=nxt
        print(f"    round {r+1}: {len(cur):,}"); sys.stdout.flush()
        if len(cur)>cap: break
    return cur

def edges(cur, units):
    E=0
    for k,(re,im) in cur.items():
        for (ure,uim) in units:
            q=(re[0]+ure[0],re[1]+ure[1],re[2]+ure[2],re[3]+ure[3],
               im[0]+uim[0],im[1]+uim[1],im[2]+uim[2],im[3]+uim[3])
            if q in cur: E+=1
    return E//2

if __name__=='__main__':
    CAP=int(sys.argv[1]) if len(sys.argv)>1 else 700_000
    UA=dirs_from(W3,2); UB=dirs_from(W4,2)
    for name,dirs,depth in (("A  Moser-angle translation",UA,7),
                            ("B  deGrey-angle translation",UB,7)):
        t=time.time(); cur=translation_closure(dirs,depth,CAP); tb=time.time()-t
        n=len(cur)
        t=time.time(); m=edges(cur,dirs); te=time.time()-t
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6
        print(f"{name}: |dirs|={len(dirs)}  n={n:,}  m={m:,}  "
              f"meandeg={2*m/n:.2f}  build {tb:.0f}s edges {te:.0f}s  RSS {rss:.2f}GB")
        sys.stdout.flush()
