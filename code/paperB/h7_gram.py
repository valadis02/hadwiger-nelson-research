"""
h7_gram.py — step D3.b3(a), the exact lift, done through the Gram matrix.

Working with coordinates wastes effort on the rotation gauge.  The Gram matrix
carries the same information with none of it:

    G is 7x7 symmetric, diagonal 1, the 11 edges force G_ij = 1/2, ten entries
    unknown, and realisability in R^3 is exactly rank(G) <= 3, i.e. every 4x4
    minor vanishes.

We refine the numerical G to 40 digits by Newton on the minor equations, then
identify each unknown entry with PSLQ.  Nine digits was not enough: PSLQ
returned a bogus quadratic with six-figure coefficients.
"""
import itertools, sys
import numpy as np
from mpmath import mp, mpf, matrix, lu_solve, pslq, nstr

mp.dps = 60

V = np.load("h7_realisation.npy")
n = 7
G0 = V @ V.T
E = [(i, j) for i in range(n) for j in range(i+1, n)
     if abs(G0[i, j] - 0.5) < 1e-9]
NE = [(i, j) for i in range(n) for j in range(i+1, n) if (i, j) not in E]
print(f"edges (inner product 1/2): {len(E)}   unknown entries: {len(NE)}")

def build(x):
    """symmetric matrix from the unknown entries"""
    M = matrix(n, n)
    for i in range(n):
        M[i, i] = mpf(1)
    for (i, j) in E:
        M[i, j] = M[j, i] = mpf(1)/2
    for t, (i, j) in enumerate(NE):
        M[i, j] = M[j, i] = x[t]
    return M

def minors(x):
    """all 4x4 principal-ish minors; rank <= 3 iff these vanish"""
    M = build(x)
    out = []
    for q in itertools.combinations(range(n), 4):
        S = matrix(4, 4)
        for a in range(4):
            for b in range(4):
                S[a, b] = M[q[a], q[b]]
        out.append(S[0,0]*(S[1,1]*(S[2,2]*S[3,3]-S[2,3]*S[3,2])
                           -S[1,2]*(S[2,1]*S[3,3]-S[2,3]*S[3,1])
                           +S[1,3]*(S[2,1]*S[3,2]-S[2,2]*S[3,1]))
                   -S[0,1]*(S[1,0]*(S[2,2]*S[3,3]-S[2,3]*S[3,2])
                           -S[1,2]*(S[2,0]*S[3,3]-S[2,3]*S[3,0])
                           +S[1,3]*(S[2,0]*S[3,2]-S[2,2]*S[3,0]))
                   +S[0,2]*(S[1,0]*(S[2,1]*S[3,3]-S[2,3]*S[3,1])
                           -S[1,1]*(S[2,0]*S[3,3]-S[2,3]*S[3,0])
                           +S[1,3]*(S[2,0]*S[3,1]-S[2,1]*S[3,0]))
                   -S[0,3]*(S[1,0]*(S[2,1]*S[3,2]-S[2,2]*S[3,1])
                           -S[1,1]*(S[2,0]*S[3,2]-S[2,2]*S[3,0])
                           +S[1,2]*(S[2,0]*S[3,1]-S[2,1]*S[3,0])))
    return out

x = [mpf(float(G0[i, j])) for (i, j) in NE]
m = len(NE)
for it in range(60):
    f = minors(x)
    F = matrix(len(f), 1)
    for i, v in enumerate(f):
        F[i] = v
    J = matrix(len(f), m)
    h = mpf(10)**(-25)
    for t in range(m):
        xp = list(x); xp[t] = xp[t] + h
        fp = minors(xp)
        for i in range(len(f)):
            J[i, t] = (fp[i] - f[i]) / h
    # least-squares Newton step via normal equations
    JT = J.T
    A = JT * J
    for i in range(m):
        A[i, i] = A[i, i] + mpf(10)**(-40)
    b = JT * F
    try:
        d = lu_solve(A, b)
    except Exception as e:
        print("  singular:", e); break
    for t in range(m):
        x[t] = x[t] - d[t]
    res = max(abs(v) for v in minors(x))
    if it % 10 == 0 or res < mpf(10)**(-45):
        print(f"  iter {it}: max |minor| = {nstr(res, 5)}")
        sys.stdout.flush()
    if res < mpf(10)**(-45):
        break

print("\nrefined non-edge inner products, and their minimal polynomials:")
seen = []
for t, (i, j) in enumerate(NE):
    v = x[t]
    if any(abs(v - w) < mpf(10)**(-30) for w in seen):
        continue
    seen.append(v)
    lin = pslq([v, mpf(1)], tol=mpf(10)**(-40), maxcoeff=10**6, maxsteps=50000)
    quad = pslq([v*v, v, mpf(1)], tol=mpf(10)**(-40), maxcoeff=10**8, maxsteps=80000)
    quart = pslq([v**4, v**3, v*v, v, mpf(1)], tol=mpf(10)**(-38),
                 maxcoeff=10**6, maxsteps=80000)
    print(f"  {nstr(v, 25)}")
    print(f"     rational: {lin}   quadratic: {quad}   quartic: {quart}")
    sys.stdout.flush()
