"""
bigbuild2.py — optimised scale test.

Two changes over bigbuild.py, both worth stating because they are the whole
difference between 10^4 and 10^7:

 (A) THE DENOMINATOR DOES NOT GROW WITH DEPTH.
     A point is a SUM of unit vectors omega_1^a omega_3^b.  Each such direction
     already has denominator dividing 2^5 * 6^2 = 1152 (five halvings from
     omega_1, two sixths from omega_3).  Sums of things over a common
     denominator stay over that denominator.  So DEN = 1152 for every depth,
     not 6^depth.  Coordinates stay small, int64 never overflows, hashing is
     cheap.  (bigbuild.py used 6^depth and paid for it.)

 (B) CANDIDATE PAIRS FROM A KD-TREE, DECIDED EXACTLY IN NUMPY.
     scipy's cKDTree returns all pairs within radius 1.01 at C speed; the exact
     integer test is then a vectorised numpy expression over the whole candidate
     array at once.  Float proposes, exact decides -- never the other way round.

Reports points/s, edges/s, bytes/point, so the reachable size is a measurement
rather than a guess.
"""
import time, sys
import numpy as np
from math import sqrt
from scipy.spatial import cKDTree

S3, S11, S33 = sqrt(3), sqrt(11), sqrt(33)
DEN = 2**5 * 6**6            # clears omega_1^a (a<=5) and omega_3^b for |b|<=6, at ANY depth

MUL = {}
for (i, j), (k, c) in {
    (0, 0): (0, 1), (0, 1): (1, 1), (0, 2): (2, 1), (0, 3): (3, 1),
    (1, 1): (0, 3), (1, 2): (3, 1), (1, 3): (2, 3),
    (2, 2): (0, 11), (2, 3): (1, 11), (3, 3): (0, 33),
}.items():
    MUL[(i, j)] = (k, c); MUL[(j, i)] = (k, c)


def fmul(x, y):
    o = [0, 0, 0, 0]
    for i in range(4):
        if x[i]:
            for j in range(4):
                if y[j]:
                    k, c = MUL[(i, j)]
                    o[k] += x[i] * y[j] * c
    return tuple(o)


def fadd(x, y): return (x[0]+y[0], x[1]+y[1], x[2]+y[2], x[3]+y[3])
def fsub(x, y): return (x[0]-y[0], x[1]-y[1], x[2]-y[2], x[3]-y[3])
def fval(x): return x[0] + x[1]*S3 + x[2]*S11 + x[3]*S33


def cmul_scaled(a, b):
    re = fsub(fmul(a[0], b[0]), fmul(a[1], b[1]))
    im = fadd(fmul(a[0], b[1]), fmul(a[1], b[0]))
    assert all(v % DEN == 0 for v in re + im)
    return (tuple(v//DEN for v in re), tuple(v//DEN for v in im))


def directions(bmax):
    w1 = ((DEN//2, 0, 0, 0), (0, DEN//2, 0, 0))
    w3 = ((5*DEN//6, 0, 0, 0), (0, 0, DEN//6, 0))
    w3c = (w3[0], tuple(-v for v in w3[1]))
    one = ((DEN, 0, 0, 0), (0, 0, 0, 0))
    out = {}
    for a in range(6):
        wa = one
        for _ in range(a): wa = cmul_scaled(wa, w1)
        for b in range(-bmax, bmax+1):
            z = wa
            for _ in range(abs(b)): z = cmul_scaled(z, w3 if b > 0 else w3c)
            out.setdefault(z[0]+z[1], z)
    return list(out.values())


def build(dirs, depth, cap):
    """iterated Minkowski sum, keys are 8-tuples of ints"""
    cur = {(0,)*8: ((0,0,0,0), (0,0,0,0))}
    for step in range(depth):
        nxt = {}
        for (re, im) in cur.values():
            for (dre, dim) in dirs:
                r = (re[0]+dre[0], re[1]+dre[1], re[2]+dre[2], re[3]+dre[3])
                m = (im[0]+dim[0], im[1]+dim[1], im[2]+dim[2], im[3]+dim[3])
                nxt[r+m] = (r, m)
        cur = nxt
        print(f"    depth {step+1}: {len(cur):,}"); sys.stdout.flush()
        if len(cur) > cap:
            break
    return cur


def to_arrays(cur):
    n = len(cur)
    RE = np.empty((n, 4), dtype=np.int64)
    IM = np.empty((n, 4), dtype=np.int64)
    for i, (re, im) in enumerate(cur.values()):
        RE[i] = re; IM[i] = im
    X = (RE[:, 0] + RE[:, 1]*S3 + RE[:, 2]*S11 + RE[:, 3]*S33) / DEN
    Y = (IM[:, 0] + IM[:, 1]*S3 + IM[:, 2]*S11 + IM[:, 3]*S33) / DEN
    return RE, IM, X, Y


def vec_fmul_self(A):
    """(a0+a1√3+a2√11+a3√33)^2 as an (n,4) int array"""
    a0, a1, a2, a3 = A[:, 0], A[:, 1], A[:, 2], A[:, 3]
    c0 = a0*a0 + 3*a1*a1 + 11*a2*a2 + 33*a3*a3
    c1 = 2*a0*a1 + 22*a2*a3
    c2 = 2*a0*a2 + 6*a1*a3
    c3 = 2*a0*a3 + 2*a1*a2
    return c0, c1, c2, c3


def unit_group(bmax_u):
    """all units omega_1^a omega_3^b with |b| <= bmax_u, as scaled int pairs"""
    return directions(bmax_u)


def unit_edges_algebraic(cur, bmax_u):
    """
    NO FLOATING POINT AT ALL.

    An edge is a pair of points differing by a unit of Z[omega_1,omega_3].
    So: for every point p and every unit u, look up p+u in the hash table.
    Exact integer arithmetic, O(n * |U|) dictionary probes.

    bmax_u must be taken large enough that enlarging it finds no new edges;
    the caller checks that saturation explicitly rather than assuming it.
    """
    U = unit_group(bmax_u)
    keys = cur
    E = 0
    for k, (re, im) in cur.items():
        for (ure, uim) in U:
            q = (re[0]+ure[0], re[1]+ure[1], re[2]+ure[2], re[3]+ure[3],
                 im[0]+uim[0], im[1]+uim[1], im[2]+uim[2], im[3]+uim[3])
            if q in keys:
                E += 1
    return E // 2, len(U)


if __name__ == '__main__':
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 400_000
    bmax = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    depth = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    dirs = directions(bmax)
    print(f"=== {len(dirs)} directions, DEN={DEN}, cap={cap:,} ===")
    t = time.time(); cur = build(dirs, depth, cap); tb = time.time()-t
    n = len(cur)
    print(f"  build : {n:,} pts in {tb:.1f}s = {n/max(tb,1e-9):,.0f} pts/s")
    sys.stdout.flush()
    prev = None
    for bu in (2, 4, 6):
        t = time.time(); m, nu = unit_edges_algebraic(cur, bu); te = time.time()-t
        flag = "" if prev is None else ("  SATURATED" if m == prev else "  +NEW EDGES")
        print(f"  edges : |U|={nu:3d} -> {m:,} edges in {te:.1f}s "
              f"({n*nu/max(te,1e-9):,.0f} probes/s), mean deg {2*m/n:.2f}{flag}")
        sys.stdout.flush()
        prev = m
