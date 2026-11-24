"""Extract a small 6-chromatic subgraph (k=5 UNSAT core) of the depth-2 ball in Q(i,sqrt3,sqrt7)."""
import sys, time, pickle, threading
import numpy as np, networkx as nx
from pysat.solvers import Solver
import ball37 as B

def unsat_core(G, k, verts, tl=120, solver="cadical195"):
    """Return a subset of verts whose induced subgraph is not k-colourable (UNSAT core), or None."""
    H = G.subgraph(verts)
    nodes = list(H.nodes()); ix = {v: i for i, v in enumerate(nodes)}; n = len(nodes)
    var = lambda v, c: ix[v] * k + c + 1
    sel = lambda v: n * k + ix[v] + 1
    s = Solver(name=solver)
    for v in nodes: s.add_clause([-sel(v)] + [var(v, c) for c in range(k)])
    for u, v in H.edges():
        for c in range(k): s.add_clause([-var(u, c), -var(v, c)])
    # symmetry breaking on a clique: only valid if the clique stays; use it only for colours of a triangle-free start -> skip
    r = s.solve(assumptions=[sel(v) for v in nodes])
    if r:
        return None
    core = s.get_core()
    s.delete()
    inv = {sel(v): v for v in nodes}
    return [inv[l] for l in core if l in inv]

def k_colourable(G, verts, k, solver="cadical195"):
    H = G.subgraph(verts)
    nodes = list(H.nodes()); ix = {v: i for i, v in enumerate(nodes)}
    var = lambda v, c: ix[v] * k + c + 1
    s = Solver(name=solver)
    for v in nodes: s.add_clause([var(v, c) for c in range(k)])
    for u, v in H.edges():
        for c in range(k): s.add_clause([-var(u, c), -var(v, c)])
    r = s.solve(); s.delete()
    return r

if __name__ == "__main__":
    U = B.units(1); gens = U + [(1 + B.zeta) * u for u in U]
    P = B.ball(gens, 2); G = B.graph(P)
    G = nx.k_core(G, 5)
    t0 = time.time()
    verts = list(G.nodes())
    for it in range(6):
        c = unsat_core(G, 5, verts)
        print(f"iter {it}: core size {len(c)} ({time.time()-t0:.1f}s)", flush=True)
        c = list(nx.k_core(G.subgraph(c), 5).nodes())
        if len(c) >= len(verts): break
        verts = c
    pickle.dump({"P": P, "verts": verts, "edges": list(G.subgraph(verts).edges())}, open("core6_stage1.pkl", "wb"))
    # deletion-based minimisation (time-boxed)
    i = 0
    order = sorted(verts, key=lambda v: G.subgraph(verts).degree(v))
    for v in order:
        if time.time() - t0 > 250: break
        if v not in verts: continue
        trial = [w for w in verts if w != v]
        trial = list(nx.k_core(G.subgraph(trial), 5).nodes())
        if not k_colourable(G, trial, 5):
            verts = trial
    H = G.subgraph(verts)
    print(f"final: |V|={H.number_of_nodes()} |E|={H.number_of_edges()} ({time.time()-t0:.1f}s)", flush=True)
    pickle.dump({"P": P, "verts": verts, "edges": list(H.edges())}, open("core6.pkl", "wb"))
