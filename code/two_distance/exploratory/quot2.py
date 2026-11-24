"""Idea 2: finite-field quotients for the two-distance problem chi({1,d}).

If p = 3 mod 4 and the prime of K+ above p is inert in K = K+(i), every finite
{1,d}-graph with vertices in K maps homomorphically to
    H(p, c) = Cay(F_p(i), {x : N(x) = 1} U {x : N(x) = c}),  c = d^2 mod P.
So chi(H) is an upper bound for every such graph.
"""
import sys, time, itertools
import networkx as nx
from pysat.solvers import Glucose4
from pysat.card import CardEnc

def H(p, c):
    V = [(a, b) for a in range(p) for b in range(p)]
    idx = {v: i for i, v in enumerate(V)}
    C = [(a, b) for (a, b) in V if (a * a + b * b) % p in ({1, c % p})]
    assert (0, 0) not in C
    G = nx.Graph(); G.add_nodes_from(range(len(V)))
    for (a, b) in V:
        for (x, y) in C:
            G.add_edge(idx[(a, b)], idx[((a + x) % p, (b + y) % p)])
    return G, len(C)

def colourable(G, k, clique, tlimit):
    n = G.number_of_nodes()
    var = lambda v, c: v * k + c + 1
    s = Glucose4()
    for v in range(n):
        s.add_clause([var(v, c) for c in range(k)])
    for u, v in G.edges():
        for c in range(k):
            s.add_clause([-var(u, c), -var(v, c)])
    for i, v in enumerate(clique[:k]):          # symmetry breaking
        s.add_clause([var(v, i)])
    import threading
    timer = threading.Timer(tlimit, lambda: s.interrupt()); timer.start()
    r = s.solve_limited(expect_interrupt=True)
    timer.cancel()
    model = s.get_model() if r else None
    s.delete()
    if r:
        col = {v: next(c for c in range(k) if model[var(v, c) - 1] > 0) for v in range(n)}
        assert all(col[u] != col[v] for u, v in G.edges())
        return True
    return False if r is False else None

def sqrt_mod(a, p):
    return [x for x in range(p) if x * x % p == a % p]

if __name__ == "__main__":
    dname = sys.argv[1]
    plist = [int(x) for x in sys.argv[2].split(",")]
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 9
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 60
    for p in plist:
        assert p % 4 == 3
        for s in sqrt_mod(3, p)[:2]:
            c = (2 + s) % p if dname == "d" else 3 % p
            G, deg = H(p, c)
            cl = max(nx.find_cliques(G), key=len) if p <= 50 else next(nx.find_cliques(G))
            res = []
            t0 = time.time()
            for k in range(len(cl), kmax + 1):
                r = colourable(G, k, cl, tl)
                res.append((k, {True: "SAT", False: "UNSAT", None: "timeout"}[r]))
                if r: break
            print(f"d={dname} p={p} sqrt3={s} c={c} |V|={p*p} deg={deg} omega>={len(cl)} {res} ({time.time()-t0:.1f}s)", flush=True)
            if dname != "d": break   # c independent of s for sqrt3
