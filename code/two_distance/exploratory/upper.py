"""Periodic tile colourings avoiding distances 1 and d.
Tiles: Voronoi cells of a lattice L (scaled), coloured by L / M, M a sublattice of index k.
Valid (with margin eps) iff diam(tile) < 1 - eps and, for every v in M \\ {0},
[dmin(T, T+v) - eps, dmax(T, T+v) + eps] contains neither 1 nor d."""
import math, itertools, sys
import numpy as np

D = (math.sqrt(2) + math.sqrt(6)) / 2
EPS = 1e-6

def voronoi_cell(B):
    """Voronoi cell of lattice with basis rows B, as polygon (counter-clockwise vertices)."""
    from scipy.spatial import Voronoi
    pts = np.array([i * B[0] + j * B[1] for i in range(-3, 4) for j in range(-3, 4)])
    vor = Voronoi(pts)
    c = np.argmin(np.linalg.norm(pts, axis=1))
    reg = vor.regions[vor.point_region[c]]
    V = vor.vertices[reg]
    ang = np.arctan2(V[:, 1], V[:, 0]); V = V[np.argsort(ang)]
    return V

def seg_dist(p, a, b):
    ab = b - a; t = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0, 1)
    return np.linalg.norm(p - (a + t * ab))

def poly_min_dist(P, Q):
    # convex, disjoint interiors assumed; min over vertex-edge pairs (0 if touching)
    m = math.inf
    for X, Y in ((P, Q), (Q, P)):
        for p in X:
            for i in range(len(Y)):
                m = min(m, seg_dist(p, Y[i], Y[(i + 1) % len(Y)]))
    return m

def poly_max_dist(P, Q):
    return max(np.linalg.norm(p - q) for p in P for q in Q)

def sublattices(k):
    for a in range(1, k + 1):
        if k % a: continue
        c = k // a
        for b in range(c):
            yield np.array([[a, b], [0, c]])   # rows: generators in lattice coordinates

def check(B, s, Mcoef, cell):
    """B lattice basis (unit scale), s scale, Mcoef 2x2 integer rows giving M generators."""
    T = cell * s
    diam = poly_max_dist(T, T)
    if diam >= 1 - EPS: return False
    Mg = Mcoef @ (B * s)
    R = 3 + D  # only vectors whose tiles can be within distance D + margin
    worst = None
    rng = int(math.ceil((R + 2 * diam) / min(np.linalg.norm(Mg[0]), np.linalg.norm(Mg[1]), 1e-9) )) + 3
    for i in range(-rng, rng + 1):
        for j in range(-rng, rng + 1):
            if i == 0 and j == 0: continue
            v = i * Mg[0] + j * Mg[1]
            if np.linalg.norm(v) - diam > D + 1: continue
            lo = poly_min_dist(T, T + v); hi = poly_max_dist(T, T + v)
            for t in (1.0, D):
                if lo - EPS <= t <= hi + EPS: return False
    return True

def search(kmax=20, steps=400):
    fams = {"hex": np.array([[1.0, 0.0], [0.5, math.sqrt(3) / 2]]), "square": np.array([[1.0, 0.0], [0.0, 1.0]])}
    for k in range(7, kmax + 1):
        print("k =", k, flush=True)
        for name, B in fams.items():
            cell = voronoi_cell(B)
            diam1 = poly_max_dist(cell, cell)
            for Mcoef in sublattices(k):
                for s in np.linspace(0.3 / diam1, (1 - 2 * EPS) / diam1, steps):
                    if check(B, s, Mcoef, cell):
                        return k, name, Mcoef, s, cell * s
    return None

if __name__ == "__main__":
    r = search(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
    print(r)
