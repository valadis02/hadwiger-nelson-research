"""Exact arithmetic in Q(sqrt3, sqrt7) and in K = Q(i, sqrt3, sqrt7) (as pairs re, im).

An element of Q(sqrt3, sqrt7) is a 4-tuple of Fractions (a, b, c, e) meaning
a + b*sqrt3 + c*sqrt7 + e*sqrt21.  No floating point is used anywhere."""
from fractions import Fraction as Fr

def R(a=0, b=0, c=0, e=0):
    return (Fr(a), Fr(b), Fr(c), Fr(e))

ZERO, ONE = R(0), R(1)

def radd(x, y): return tuple(p + q for p, q in zip(x, y))
def rsub(x, y): return tuple(p - q for p, q in zip(x, y))
def rneg(x): return tuple(-p for p in x)

def rmul(x, y):
    a, b, c, e = x; A, B, C, E = y
    # basis 1, s3, s7, s21 with s3^2=3, s7^2=7, s21^2=21, s3*s7=s21, s3*s21=3*s7, s7*s21=7*s3
    return (a*A + 3*b*B + 7*c*C + 21*e*E,
            a*B + b*A + 7*(c*E + e*C),
            a*C + c*A + 3*(b*E + e*B),
            a*E + e*A + b*C + c*B)

def rscale(x, q): return tuple(p * Fr(q) for p in x)

# complex numbers over Q(sqrt3, sqrt7)
def C(re, im=ZERO): return (re, im)
def cadd(z, w): return (radd(z[0], w[0]), radd(z[1], w[1]))
def csub(z, w): return (rsub(z[0], w[0]), rsub(z[1], w[1]))
def cmul(z, w): return (rsub(rmul(z[0], w[0]), rmul(z[1], w[1])), radd(rmul(z[0], w[1]), rmul(z[1], w[0])))
def cconj(z): return (z[0], rneg(z[1]))
def norm2(z): return radd(rmul(z[0], z[0]), rmul(z[1], z[1]))   # |z|^2, an element of Q(sqrt3, sqrt7)

CZERO, CONE = C(ZERO), C(ONE)
ZETA = C(R(0, Fr(1, 2)), R(Fr(1, 2)))            # e^{i pi/6} = (sqrt3 + i)/2
OMEGA = C(R(Fr(3, 4)), R(0, 0, Fr(1, 4)))        # (3 + i sqrt7)/4

def cpow(z, k):
    r = CONE
    for _ in range(k): r = cmul(r, z)
    return r

D2 = R(2, 1)          # d^2 = 2 + sqrt3,  d = (sqrt2 + sqrt6)/2
FORBIDDEN = (ONE, D2)

def generators():
    """index 3k + (a+1) = zeta^k omega^a (k = 0..11, a = -1, 0, 1); index 36 + m = (1 + zeta) * gen[m]."""
    g = []
    for k in range(12):
        for a in (-1, 0, 1):
            w = cconj(OMEGA) if a == -1 else (CONE if a == 0 else OMEGA)
            g.append(cmul(cpow(ZETA, k), w))
    one_zeta = cadd(CONE, ZETA)
    g += [cmul(one_zeta, u) for u in g[:36]]
    return g

def selfcheck():
    assert norm2(ZETA) == ONE and norm2(OMEGA) == ONE
    assert cpow(ZETA, 12) == CONE
    assert norm2(cadd(CONE, ZETA)) == D2
    s3 = R(0, 1); assert rmul(s3, s3) == R(3)
    s7 = R(0, 0, 1); assert rmul(s7, s7) == R(7)
    assert rmul(s3, s7) == R(0, 0, 0, 1)

selfcheck()
