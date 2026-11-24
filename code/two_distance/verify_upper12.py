"""chi({1,d}) <= 12: exact check of the separation of the hexagonal 12-colouring.

Tile centres: triangular lattice with unit spacing (e1 = 1, e2 = e^{i pi/3}); tiles are the Voronoi
hexagons, of diameter 2/sqrt3.  Colour classes are cosets of M = <2(e1+e2), 6 e2> (index 12).
Two tiles T, T+v of the same colour are at distance >= |v| - diam.  We check exactly that
min_{v in M, v != 0} |v|^2 = 12, so the separation is >= 2*sqrt3 - 2/sqrt3 = 2 * diam.
After scaling so that diam = D with d/2 < D < 1, every distance in [D, 2D) - in particular 1 and d -
is avoided."""
from fractions import Fraction as Fr
import itertools
# |x e1 + y e2|^2 = x^2 + x y + y^2 for the triangular lattice
q = lambda x, y: x * x + x * y + y * y
gens = [(2, 2), (0, 6)]
det = abs(gens[0][0] * gens[1][1] - gens[0][1] * gens[1][0])
assert det == 12
best = min(q(i * 2 + j * 0, i * 2 + j * 6) for i, j in itertools.product(range(-6, 7), repeat=2) if (i, j) != (0, 0))
assert best == 12, best
# 2*sqrt3 - 2/sqrt3 = 4/sqrt3 = 2 * (2/sqrt3): squared, 16/3 = 4 * 4/3
assert Fr(16, 3) == 4 * Fr(4, 3)
print("index 12, min |v|^2 = 12, separation = 2 x tile diameter  OK")
d = (2 ** 0.5 + 6 ** 0.5) / 2
assert d / 2 < 0.98 < 1 and 1 < 2 * 0.98 and d < 2 * 0.98
print(f"with D = 0.98: avoided interval [0.98, 1.96) contains 1 and d = {d:.4f}  OK")
