"""
build_g17.py — ancillary file for "Monochromatic-pair probabilities for
5-colourings of the plane", section 4.2.

Builds the Parshall / Mixon / Alexeev graph in exact arithmetic over
Q(sqrt3, sqrt11) and asserts every count stated in the text:
    19 vertices / 45 unit edges   (the three-hexagon union)
    21 vertices                   (Parshall)
    17 vertices / 38 unit edges   (after Mixon's four deletions)
    98 non-unit pairs in 25 distance classes, 15 at sqrt3 and 7 at 2,
    and no two classes numerically equal.

Usage: python3 build_g17.py
"""
import itertools
import sympy as sp

I = sp.I
w1 = (1 + I*sp.sqrt(3))/2                      # e^{i pi/3}
w3 = (sp.Integer(5) + I*sp.sqrt(11))/6         # the Moser angle, cos = 5/6
hexs = [sp.Integer(0), sp.Integer(1), w1, 1 + w1,
        sp.conjugate(w1), 1 + sp.conjugate(w1), sp.Integer(2)]
v = 1 + w1


def dedupe(points):
    out = []
    for p in points:
        p = sp.expand(sp.simplify(p))
        if not any(sp.simplify(p - q) == 0 for q in out):
            out.append(p)
    return out


def unit_edges(points):
    return [(a, b) for a, b in itertools.combinations(range(len(points)), 2)
            if sp.simplify(sp.expand(sp.Abs(points[a] - points[b])**2)) == 1]


# ---- 19 vertices / 45 edges --------------------------------------------
green = [w3*p for p in hexs]
blue = [v + w3*(p - v) for p in hexs]
G19 = dedupe(hexs + green + blue)
assert len(G19) == 19, len(G19)
assert len(unit_edges(G19)) == 45, len(unit_edges(G19))
print("three-hexagon union : 19 vertices, 45 unit edges   OK")

# ---- Parshall's 21 ------------------------------------------------------
G21 = dedupe(hexs + green + blue + [2 + w1, sp.Integer(3)])
assert len(G21) == 21, len(G21)
print("Parshall            : 21 vertices                  OK")

# ---- Mixon's four deletions -> 17 --------------------------------------
rem = [sp.Integer(3),
       sp.Rational(1, 2)*I*sp.sqrt(3) + sp.Rational(5, 2),
       sp.Rational(1, 12)*(I*sp.sqrt(11) - 1)*(-I*sp.sqrt(3) - 3),
       sp.Rational(1, 12)*(-I*sp.sqrt(11) - 5)*(-I*sp.sqrt(3) - 1)
       + sp.Rational(1, 12)*(I*sp.sqrt(11) - 1)*(-I*sp.sqrt(3) - 3)]
rem = [sp.expand(sp.simplify(r)) for r in rem]
assert all(any(sp.simplify(r - q) == 0 for q in G21) for r in rem), \
    "each deleted point must be a vertex of the 21-vertex graph"

pts = [p for p in G21 if not any(sp.simplify(p - r) == 0 for r in rem)]
n = len(pts)
assert n == 17, n
E = unit_edges(pts)
assert len(E) == 38, len(E)
print("after deletions     : 17 vertices, 38 unit edges   OK")

# ---- distance classes ---------------------------------------------------
d2 = {}
for a, b in itertools.combinations(range(n), 2):
    d2[(a, b)] = sp.simplify(sp.expand(sp.Abs(pts[a] - pts[b])**2))

classes = {}
for e, d in d2.items():
    if d != 1:
        classes.setdefault(sp.srepr(d), []).append(e)

names = sorted(classes)
assert sum(len(c) for c in classes.values()) == 98
assert len(names) == 25, len(names)
assert len(classes[sp.srepr(sp.Integer(3))]) == 15      # distance sqrt3
assert len(classes[sp.srepr(sp.Integer(4))]) == 7       # distance 2

vals = sorted((float(sp.N(sp.sympify(k))), k) for k in names)
assert all(abs(vals[i][0] - vals[i + 1][0]) > 1e-12 for i in range(len(vals) - 1)), \
    "two distance classes coincide numerically - symbolic forms are not canonical"
print("distance classes    : 98 pairs in 25 classes       OK")
print("                      15 at sqrt3, 7 at 2, no collisions")

print("\n  d^2                          pairs")
for val, nm in vals:
    print(f"  {str(sp.sympify(nm)):26s}  {len(classes[nm]):3d}")

print("\nall assertions passed")
