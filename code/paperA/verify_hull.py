"""
verify_hull.py — ancillary file for "Monochromatic-pair probabilities for
5-colourings of the plane", section 4.1.

Enumerates all proper k-colourings of H = W_6 for k = 4, 5; computes the exact
points (a/6, b/3), their Pareto minima and the facets of the lower-left hull;
derives Theorems 4.1-4.4 in rational arithmetic and asserts every stated constant.

Fractions only. No floating point anywhere.

Usage: python3 verify_hull.py
"""
from fractions import Fraction as F
from itertools import product

# ---- H = W_6: centre 0, rim 1..6 in a cycle -----------------------------
rim = list(range(1, 7))
edges = [(0, i) for i in rim] + [(i, 1 + (i % 6)) for i in rim]
sqrt3_pairs = sorted({tuple(sorted((i, 1 + ((i + 1) % 6)))) for i in rim})
two_pairs = sorted({tuple(sorted((i, 1 + ((i + 2) % 6)))) for i in rim})
assert len(edges) == 12
assert len(sqrt3_pairs) == 6 and len(two_pairs) == 3
assert all(0 not in p for p in sqrt3_pairs + two_pairs), "centre must not occur"


def states(k):
    """exact points (a/6, b/3) over all proper k-colourings of H"""
    pts = set()
    for col in product(range(k), repeat=7):
        if any(col[u] == col[v] for u, v in edges):
            continue
        a = sum(1 for u, v in sqrt3_pairs if col[u] == col[v])
        b = sum(1 for u, v in two_pairs if col[u] == col[v])
        pts.add((F(a, 6), F(b, 3)))
    return sorted(pts)


def hull(pts):
    """Pareto minima and facets alpha*x + beta*y >= gamma of the lower-left hull"""
    P = sorted(set(pts))
    par = sorted(p for p in P
                 if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in P))
    facets = []
    for (x1, y1), (x2, y2) in zip(par, par[1:]):
        al, be = y1 - y2, x2 - x1              # outward normal, both >= 0
        if al < 0 or be < 0:
            continue
        ga = al * x1 + be * y1
        if not all(al * x + be * y >= ga for x, y in P):
            continue
        g = min(z for z in (al, be) if z)      # normalise
        f = (al / g, be / g, ga / g)
        if f not in facets:
            facets.append(f)
    return par, facets


report = {}
for k in (4, 5):
    pts = states(k)
    par, fac = hull(pts)
    report[k] = (pts, par, fac)
    print(f"\n=== H, k = {k} ===")
    print(f"  distinct (a/6, b/3) : {len(pts)}")
    print(f"  Pareto minima       : {[(str(x), str(y)) for x, y in par]}")
    for al, be, ga in fac:
        print(f"  facet               : {al}*p_sqrt3 + {be}*p_2 >= {ga}")

# ---- the table of section 4.1 ------------------------------------------
pts4, par4, fac4 = report[4]
pts5, par5, fac5 = report[5]

assert len(pts4) == 4, "k=4 must give the four states of de Grey's Figure 1"
assert len(pts5) == 8
assert par5 == [(F(0), F(2, 3)), (F(1, 6), F(1, 3)), (F(1, 3), F(0))]
(x1, y1), (x2, y2), (x3, y3) = par5
assert (x2 - x1) * (y3 - y1) == (x3 - x1) * (y2 - y1), "k=5 frontier must be collinear"
assert len(fac4) == 2 and len(fac5) == 1, "two facets at k=4, one at k=5"
assert (F(2), F(1), F(1)) in fac4          # 2 p_sqrt3 + p_2 >= 1   (c_4 = 1)
assert (F(1), F(1), F(2, 3)) in fac4       # p_sqrt3 + p_2 >= 2/3   (Lemma 15)

# ---- Theorem 4.1 -------------------------------------------------------
A, B, G = fac5[0]
assert (A, B, G) == (F(2), F(1), F(2, 3))
print(f"\nTheorem 4.1 : {A} p_sqrt3 + {B} p_2 >= {G}          [unique facet]")

assert F(1, 6) + F(1, 3) < F(2, 3), "state (1/6,1/3) must violate p_sqrt3 + p_2 >= 2/3"
print("            and p_sqrt3 + p_2 >= 2/3 is FALSE at k=5, witness (1/6, 1/3)")

# ---- Theorem 4.2 : Lemma 2 gives p_2 <= 1/2 ----------------------------
p2_max = F(1, 2)
thm42 = (G - B * p2_max) / A
assert thm42 == F(1, 12)
print(f"Theorem 4.2 : p^(5)_sqrt3 >= {thm42}")

# ---- Theorem 4.3 : Lemma 38 on (3,4,7) and (3,7,7) ---------------------
#   p_sqrt3 + p_2 + p_sqrt7 <= 1
#   p_sqrt3        + 2 p_sqrt7 <= 1
#   sum: 2 p_sqrt3 + p_2 + 3 p_sqrt7 <= 2 ; with Theorem 4.1, 3 p_sqrt7 <= 2 - 2/3
thm43 = (F(2) - G) / 3
assert thm43 == F(4, 9)
print(f"Theorem 4.3 : p^(5)_sqrt7 <= {thm43}")

# ---- Theorem 4.4 : same, at k=4, with Lemma 15 -------------------------
lemma15 = F(2, 3)
thm44 = F(1) - lemma15
assert thm44 == F(1, 3)
assert thm44 < F(37, 100), "must improve the wiki's 37/100"
print(f"Theorem 4.4 : p^(4)_sqrt7 <= {thm44}   (wiki: 37/100)")

print("\nall assertions passed")
