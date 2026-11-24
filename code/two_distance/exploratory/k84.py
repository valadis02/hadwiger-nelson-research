"""Search for K8-4 gadgets (8 points, all pairs at distance 1 or d except a perfect matching
of 4 pairs, each at distance > 1/2) for a given d, in Minkowski balls of chosen rotations.
A K8-4 gadget with k=6 forces sum of the 4 matched-pair probabilities >= 2, i.e. all = 1/2."""
import sys, math, cmath, itertools, time
import numpy as np
from scipy.spatial import cKDTree
import networkx as nx

def ball(gens, depth):
    pts = {(0.0, 0.0): 0j}; fr = [0j]
    for _ in range(depth):
        new = []
        for p in fr:
            for g in gens:
                q = p + g; key = (round(q.real, 7), round(q.imag, 7))
                if key not in pts: pts[key] = q; new.append(q)
        fr = new
    return np.array(list(pts.values()))

def search(P, d2, maxhits=20):
    xy = np.c_[P.real, P.imag]; T = cKDTree(xy); n = len(P)
    G = nx.Graph(); G.add_nodes_from(range(n))
    for r2 in (1.0, d2):
        for i, j in T.query_pairs(math.sqrt(r2) + 1e-6):
            if abs(((xy[i] - xy[j]) ** 2).sum() - r2) < 1e-7: G.add_edge(i, j)
    zero = int(np.argmin(np.abs(P)))
    # pairs (0, b): b non-adjacent to 0, |b| > 1/2, with large common neighbourhood
    hits = []
    N0 = set(G[zero])
    for b in range(n):
        if b == zero or G.has_edge(zero, b) or abs(P[b]) <= 0.5: continue
        C = N0 & set(G[b])
        if len(C) < 6: continue
        # inside C find 3 disjoint non-adjacent pairs, all cross pairs adjacent: a 3-partite K_{2,2,2}
        Cl = list(C)
        pairs = [(u, v) for u, v in itertools.combinations(Cl, 2) if not G.has_edge(u, v) and abs(P[u] - P[v]) > 0.5]
        H = nx.Graph()
        for a in range(len(pairs)):
            u, v = pairs[a]
            for c in range(a + 1, len(pairs)):
                x, y = pairs[c]
                if len({u, v, x, y}) == 4 and all(G.has_edge(s, t) for s in (u, v) for t in (x, y)):
                    H.add_edge(a, c)
        for clq in nx.find_cliques(H):
            if len(clq) >= 3:
                for tri in itertools.combinations(clq, 3):
                    if all(H.has_edge(s, t) for s, t in itertools.combinations(tri, 2)):
                        X = [zero, b] + [w for a in tri for w in pairs[a]]
                        dists = sorted(round(abs(P[s] - P[t]), 6) for s, t in [(zero, b)] + [pairs[a] for a in tri])
                        hits.append((X, dists))
                        if len(hits) >= maxhits: return hits, G
    return hits, G

if __name__ == "__main__":
    dname = sys.argv[1]; depth = int(sys.argv[2]); rots = sys.argv[3]
    d2 = {"sqrt3": 3.0, "2": 4.0, "phi": ((1 + 5 ** .5) / 2) ** 2, "dd": 2 + 3 ** .5}[dname]
    d = math.sqrt(d2)
    z = cmath.exp(1j * math.pi / 6)
    R = {"z12": [z ** k for k in range(12)],
         "z12p": [z ** k * r for k in range(12) for r in (1, (3 + 4j) / 5, (3 - 4j) / 5)],
         "z12w": [z ** k * r for k in range(12) for r in (1, (3 + 1j * 7 ** .5) / 4, (3 - 1j * 7 ** .5) / 4)],
         "z12q": [z ** k * r for k in range(12) for r in (1, (1 + 1j * 15 ** .5) / 4, (1 - 1j * 15 ** .5) / 4)],
         "z10": [cmath.exp(1j * math.pi / 5 * k) for k in range(10)]}[rots]
    gens = R + [d * u for u in R]
    t0 = time.time()
    P = ball(gens, depth)
    hits, G = search(P, d2)
    print(f"d={dname} rots={rots} depth={depth} |P|={len(P)} |E|={G.number_of_edges()} K8-4 gadgets found: {len(hits)} ({time.time()-t0:.0f}s)")
    for X, ds in hits[:5]: print("   matched distances:", ds)
