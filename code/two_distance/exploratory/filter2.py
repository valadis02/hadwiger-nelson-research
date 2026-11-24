"""Idea 2 filter: for K+ = Q(sqrt3, sqrt m), K = K+(i), list the non-split primes
P of K+ (inert or ramified in K/K+) with small norm and the quotient graph that
bounds chi of every finite {1,d}-graph with vertices in K."""
import sys, itertools, time, threading
import cypari2
import networkx as nx
from pysat.solvers import Glucose4
pari = cypari2.Pari(); pari.allocatemem(2*10**9)

def sat_colourable(G, k, tl=30):
    n = G.number_of_nodes(); nodes = list(G.nodes()); ix = {v: i for i, v in enumerate(nodes)}
    var = lambda v, c: ix[v] * k + c + 1
    s = Glucose4()
    for v in nodes: s.add_clause([var(v, c) for c in range(k)])
    for u, v in G.edges():
        for c in range(k): s.add_clause([-var(u, c), -var(v, c)])
    cl = max(nx.find_cliques(G), key=len) if n <= 300 else next(nx.find_cliques(G))
    for i, v in enumerate(cl[:k]): s.add_clause([var(v, i)])
    t = threading.Timer(tl, s.interrupt); t.start()
    r = s.solve_limited(expect_interrupt=True); t.cancel()
    if r:
        m = set(l for l in s.get_model() if l > 0)
        col = {v: next(c for c in range(k) if var(v, c) in m) for v in nodes}
        assert all(col[u] != col[v] for u, v in G.edges()), "bad colouring"
    s.delete()
    return r, len(cl)

def chi_range(G, kmax=8, tl=15):
    om = len(max(nx.find_cliques(G), key=len)) if G.number_of_nodes() <= 300 else 2
    lo = om
    for k in range(om, kmax + 1):
        r, _ = sat_colourable(G, k, tl)
        if r is True: return (lo, k)
        if r is False: lo = k + 1
    return (lo, None)

CACHE = {}
def inert_quotient(p, f, cmin):
    """Cay(F_{q^2}, N^-1(1) U N^-1(c)), q = p^f; cmin = minpoly of c over F_p (pari)."""
    key = ("i", p, f, str(cmin))
    if key in CACHE: return CACHE[key]
    q = p ** f
    g = pari.ffgen(pari(p) ** (2 * f), 'a')
    els = [pari(0) * g]
    # enumerate F_{q^2} as polynomials in g of degree < 2f
    els = []
    for coeffs in itertools.product(range(p), repeat=2 * f):
        x = sum(c * g ** i for i, c in enumerate(coeffs)) + 0 * g
        els.append(x)
    roots = pari.polrootsmod(cmin, p) if False else None
    # embed c: any root of cmin in F_{q^2}
    cl = pari.lift(cmin)
    cr = next(x for x in els if pari.subst(cl, 'x', x) == 0)
    norms = [x ** (q + 1) for x in els]
    S = [i for i, x in enumerate(els) if x != 0 and (norms[i] == 1 or norms[i] == cr)]
    index = {str(x): i for i, x in enumerate(els)}
    G = nx.Graph(); G.add_nodes_from(range(len(els)))
    Sel = [els[i] for i in S]
    for i, x in enumerate(els):
        for y in Sel:
            G.add_edge(i, index[str(x + y)])
    res = chi_range(G)
    CACHE[key] = (res, len(els), len(S))
    return CACHE[key]

def analyse(m, d2, pmax=30):
    pol = pari('x^2-3') if m == 1 else pari.polcompositum(pari('x^2-3'), pari(f'x^2-{m}'))[0]
    nfp = pari.nfinit(pari.subst(pol, 'x', 'y'))
    # d2 as element: d2 = a + b*sqrt3 ; find sqrt3 in nfp
    s3 = pari.nfroots(nfp, pari('x^2-3'))[0]
    a, b = d2
    c_el = pari.nfeltadd(nfp, a, pari.nfeltmul(nfp, b, s3))
    rnf = pari.rnfinit(nfp, pari('x^2+1'))
    out = []
    for p in pari.primes(pari.primepi(pmax)):
        for P in pari.idealprimedec(nfp, p):
            above = pari.rnfidealprimedec(rnf, P)
            if len(above) != 1: continue           # splits in K
            f = int(P[3]); e_rel = int(above[0][2]) // int(P[2]) if False else None
            # ramified in K/K+ ?  compare with relative discriminant
            ram = pari.idealval(nfp, pari.rnfdisc(nfp, pari('x^2+1'))[0], P) > 0
            if pari.nfeltval(nfp, c_el, P) != 0:
                out.append((int(p), f, "ram" if ram else "inert", "c=0 unusable")); continue
            modpr = pari.nfmodprinit(nfp, P)
            cbar = pari.nfmodpr(nfp, c_el, modpr)
            cmin = pari.minpoly(cbar)
            if ram:
                if f > 1:
                    out.append((int(p), f, "ram", "not computed")); continue
                cc = int(pari.lift(pari.polcoef(pari.lift(cmin), 0))) * -1 % int(p)
                S = [x for x in range(1, int(p)) if (x*x) % p in (1, cc)]
                G = nx.Graph(); G.add_nodes_from(range(int(p)))
                for x in range(int(p)):
                    for y in S: G.add_edge(x, (x + y) % int(p))
                G.remove_edges_from(nx.selfloop_edges(G))
                if any(y % p == 0 for y in S): out.append((int(p), f, "ram", "loop!")); continue
                lo, hi = chi_range(G)
                out.append((int(p), f, "ram", f"q={p} S={S} chi in [{lo},{hi}]")); continue
            if p ** f > 16:
                out.append((int(p), f, "inert", f"q={p**f} too big")); continue
            (lo, hi), nV, deg = inert_quotient(int(p), f, cmin)
            out.append((int(p), f, "inert", f"q={p**f} |V|={nV} deg={deg} chi in [{lo},{hi}]"))
    return out

if __name__ == "__main__":
    d = {"d": (2, 1), "sqrt3": (3, 0)}[sys.argv[1]]
    for m in [int(x) for x in sys.argv[2].split(",")]:
        t0 = time.time()
        print(f"--- K+ = Q(sqrt3, sqrt{m})", flush=True)
        for row in analyse(m, d):
            print("   ", row, flush=True)
        print(f"   ({time.time()-t0:.1f}s)", flush=True)
