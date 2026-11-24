"""Joint-distribution LP for chi({1,d}) > k on a finite point set F.

If c is a proper k-colouring of the plane, the invariant mean pushes forward to a probability
distribution lambda on proper k-colourings of F (proper for distances 1, d) such that
P(phi(a)=phi(b)) is the same for all pairs {a,b} of F at the same distance.
Infeasibility of this LP proves chi({1,d}) > k.  Solved by column generation (phase-1 LP,
pricing = MILP over proper colourings)."""
import sys, math, time, itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix

def classes_of(P, forb2, tol=1e-7):
    n = len(P); E = []; cls = {}
    for i in range(n):
        for j in range(i + 1, n):
            r2 = abs(P[i] - P[j]) ** 2
            if any(abs(r2 - f) < tol for f in forb2): E.append((i, j)); continue
            key = round(r2, 6)
            cls.setdefault(key, []).append((i, j))
    return E, {k: v for k, v in cls.items() if len(v) >= 2}

def price(n, k, E, w, time_limit=60):
    """max sum_e w_e * mono_e(phi) over proper k-colourings phi. w: dict pair->weight."""
    pairs = [e for e, x in w.items() if abs(x) > 1e-12]
    nx_ = n * k; nm = len(pairs) * k; N = nx_ + nm
    X = lambda v, c: v * k + c
    M = lambda t, c: nx_ + t * k + c
    rows, lo, hi = [], [], []
    A = lil_matrix((n + len(E) * k + 2 * nm + nm, N))
    r = 0
    for v in range(n):
        for c in range(k): A[r, X(v, c)] = 1
        lo.append(1); hi.append(1); r += 1
    for u, v in E:
        for c in range(k):
            A[r, X(u, c)] = 1; A[r, X(v, c)] = 1; lo.append(-np.inf); hi.append(1); r += 1
    for t, (u, v) in enumerate(pairs):
        for c in range(k):
            A[r, M(t, c)] = 1; A[r, X(u, c)] = -1; lo.append(-np.inf); hi.append(0); r += 1
            A[r, M(t, c)] = 1; A[r, X(v, c)] = -1; lo.append(-np.inf); hi.append(0); r += 1
            A[r, M(t, c)] = 1; A[r, X(u, c)] = -1; A[r, X(v, c)] = -1; lo.append(-1); hi.append(np.inf); r += 1
    A = csr_matrix(A[:r])
    obj = np.zeros(N)
    for t, e in enumerate(pairs):
        for c in range(k): obj[M(t, c)] = -w[e]
    # symmetry: vertex 0 colour 0
    lb = np.zeros(N); ub = np.ones(N); lb[X(0, 0)] = 1
    res = milp(obj, constraints=LinearConstraint(A, lo, hi), integrality=np.ones(N), bounds=Bounds(lb, ub),
               options={"time_limit": time_limit, "disp": False})
    if res.x is None or res.status != 0: return None, None   # only proven-optimal pricing counts
    x = np.round(res.x[:nx_]).reshape(n, k)
    col = x.argmax(1)
    assert all(col[u] != col[v] for u, v in E)
    return col, -res.fun

def mono_vec(col, rowpairs):
    return np.array([(1.0 if col[a] == col[b] else 0.0) - (1.0 if col[c] == col[d] else 0.0) for (a, b), (c, d) in rowpairs])

def run(P, forb2, k, iters=200, verbose=True, seed_cols=()):
    n = len(P); E, C = classes_of(P, forb2)
    rowpairs = []
    for key, prs in C.items():
        base = prs[0]
        for e in prs[1:]: rowpairs.append((e, base))
    R = len(rowpairs)
    if verbose: print(f"|F|={n} forbidden={len(E)} classes={len(C)} rows={R}", flush=True)
    cols = [np.array(c) for c in seed_cols]
    if not cols:
        c0, _ = price(n, k, E, {rowpairs[0][0]: 0.0} if R else {}, 30) if False else (None, None)
    # initial column: any proper colouring
    if not cols:
        col, _ = price(n, k, E, {}, 60)
        if col is None: print("no proper colouring at all"); return "graph not k-colourable"
        cols.append(col)
    t0 = time.time()
    for it in range(iters):
        A = np.array([mono_vec(c, rowpairs) for c in cols]).T if R else np.zeros((0, len(cols)))
        m = len(cols)
        # phase 1: min sum(s+ + s-)  s.t. A lam + s+ - s- = 0, sum lam = 1
        Aeq = np.hstack([A, np.eye(R), -np.eye(R)]); Aeq = np.vstack([Aeq, np.r_[np.ones(m), np.zeros(2 * R)]])
        beq = np.r_[np.zeros(R), 1.0]
        cost = np.r_[np.zeros(m), np.ones(2 * R)]
        res = linprog(cost, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * (m + 2 * R), method="highs")
        obj = res.fun; y = res.eqlin.marginals
        yr, y0 = y[:R], y[R]
        # reduced cost of new column a: 0 - (yr.a + y0). Want max yr.a + y0 > 0
        w = {}
        for (e1, e2), val in zip(rowpairs, yr):
            if abs(val) > 1e-12:
                w[e1] = w.get(e1, 0) + val; w[e2] = w.get(e2, 0) - val
        col, best = price(n, k, E, w, 120)
        gain = (best + y0) if best is not None else None
        if verbose: print(f"it {it}: cols={m} phase1={obj:.6f} pricing gain={gain if gain is None else round(gain,6)} ({time.time()-t0:.0f}s)", flush=True)
        if obj < 1e-9: return "FEASIBLE (no contradiction)"
        if gain is None: return "pricing failed"
        if gain <= 1e-9:
            return f"INFEASIBLE: phase-1 optimum {obj:.6f} > 0 with no improving column => chi > {k} (needs exact re-check)"
        cols.append(col)
    return "iteration limit"
