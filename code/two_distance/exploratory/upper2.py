"""Fast version: distances scale with the tile. With x = 1/s (unit distance in tile units),
need x > diam, and neither x nor D*x in any [lo_v, hi_v] (v in M\\{0}), with margin."""
import math, sys
import numpy as np
from upper import voronoi_cell, poly_min_dist, poly_max_dist, sublattices, D

def gaps(intervals, start, stop):
    """complement of union of intervals within (start, stop)."""
    iv = sorted((a, b) for a, b in intervals if b > start and a < stop)
    out = []; cur = start
    for a, b in iv:
        if a > cur: out.append((cur, a))
        cur = max(cur, b)
    if cur < stop: out.append((cur, stop))
    return out

def feasible_x(B, Mcoef, cell, XMAX=8.0, margin=1e-6):
    diam = poly_max_dist(cell, cell)
    Mg = Mcoef @ B
    lim = D * XMAX + diam
    nmin = min(np.linalg.norm(Mg[0]), np.linalg.norm(Mg[1]))
    # enumerate v with |v| <= lim + diam via a generous box in M coordinates
    G = Mg @ Mg.T; Ginv = np.linalg.inv(G)
    rng = int(math.ceil((lim + diam) * math.sqrt(max(np.linalg.eigvalsh(Ginv))))) + 1
    iv = []
    for i in range(-rng, rng + 1):
        for j in range(-rng, rng + 1):
            if i == 0 and j == 0: continue
            v = i * Mg[0] + j * Mg[1]
            if np.linalg.norm(v) - diam > lim: continue
            lo = poly_min_dist(cell, cell + v); hi = poly_max_dist(cell, cell + v)
            iv.append((lo - margin, hi + margin))
    # x must avoid iv, D*x must avoid iv  => x avoids iv and iv/D
    allv = iv + [(a / D, b / D) for a, b in iv]
    return [(a, b) for a, b in gaps(allv, diam + margin, XMAX) if b - a > 1e-9]

if __name__ == "__main__":
    fams = {"hex": np.array([[1.0, 0.0], [0.5, math.sqrt(3) / 2]]),
            "square": np.array([[1.0, 0.0], [0.0, 1.0]])}
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    found = False
    for k in range(7, kmax + 1):
        for name, B in fams.items():
            cell = voronoi_cell(B)
            for Mcoef in sublattices(k):
                g = feasible_x(B, Mcoef, cell)
                if g:
                    print(f"k={k} {name} M={Mcoef.tolist()} feasible 1/s intervals: {[(round(a,4), round(b,4)) for a,b in g]}", flush=True)
                    found = True
        if found: break
        print("k", k, "none", flush=True)
