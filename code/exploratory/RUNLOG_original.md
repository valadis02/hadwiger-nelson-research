# Run log

```
$ python3 verify_hull.py
======================================================================
CONTROLS
======================================================================
  C1  proper 4-colourings of H, mod colour permutation : 4 (expect 4)   PASS
  C2  of those, with a monochromatic sqrt3-triad       : 2 (expect 2)   PASS

======================================================================
HULLS
======================================================================
  k = 4
    colourings mod colour permutation: 11
    states (also mod symmetry)       : 4
    distinct profiles (a/6, b/3)    : 4  [('0', '1'), ('1/3', '1/3'), ('2/3', '0'), ('1', '0')]
    Pareto minima                   : [('0', '1'), ('1/3', '1/3'), ('2/3', '0')]
    facet: 2 p_sqrt3 + 1 p_2 >= 1
    facet: 3 p_sqrt3 + 3 p_2 >= 2

  k = 5
    colourings mod colour permutation: 31
    states (also mod symmetry)       : 9
    distinct profiles (a/6, b/3)    : 8  [('0', '2/3'), ('0', '1'), ('1/6', '1/3'), ('1/3', '0'), ('1/3', '1/3'), ('1/2', '0'), ('2/3', '0'), ('1', '0')]
    Pareto minima                   : [('0', '2/3'), ('1/6', '1/3'), ('1/3', '0')]
    facet: 6 p_sqrt3 + 3 p_2 >= 2

======================================================================
THEOREMS
======================================================================
  C3  k=4 hull reproduces Lemma 15 (p+p>=2/3 and 2p+p>=1) : PASS
  Thm  k=5 has a UNIQUE facet, 2 p_sqrt3 + p_2 >= 2/3      : PASS
  Thm  the three k=5 Pareto minima are COLLINEAR           : PASS
  Thm  p^(5)_sqrt3 >= 1/12   (expect 1/12)                   : PASS
  Thm  p^(5)_sqrt7 <= 4/9   (expect 4/9)                    : PASS
  Thm  p^(4)_sqrt7 <= 1/3   (expect 1/3, replacing 37/100)  : PASS

All arithmetic above is exact (fractions.Fraction); no floating
point value is computed or compared anywhere in this file.

$ python3 h7_theorem.py
========================================================================
The extremal Gram matrix over Q(sqrt46)
========================================================================
  all diagonal entries are 1 (unit vectors)      : PASS
  the 11 prescribed entries equal 1/2             : PASS
  NO non-edge entry equals 1/2                   : PASS
  rank G = 3  (expect 3)                         : PASS
  characteristic polynomial factors as:
     x**4*(32*x - 57)*(288*x**2 - 1503*x + 1790)/9216
  multiplicity of the eigenvalue 0 = 4  (expect 4)  : PASS
  57/32 and 167/64 +- sqrt(21881)/192 are roots  : PASS
  all nonzero eigenvalues >= 0  (PSD)            : PASS   [1.78125, 3.379804, 1.838946]

========================================================================
The graph realised
========================================================================
  n = 7, m = 11  (expect 11)                      : PASS
  degree sequence [3, 3, 3, 3, 3, 3, 4]  (expect [3,3,3,3,3,3,4])
  independence number = 2 (expect 2)             : PASS
  clique number       = 3 (expect 3)             : PASS
  isomorphic to the Moser spindle MS(2)          : YES

========================================================================
Cross-check: Bachoc et al. arXiv:2504.11086, Theorem 5.2
========================================================================
  polynomial 8k^2 t^3 - (k^2-10k+1) t^2 - 2(k-1) t - 1, k = 2
  roots: [-0.455258, -0.26883, 0.255338]
  t_{2,3} = 0.255338 ;  1/2 >= t_{2,3} : YES -> MS(2) is (3,1/2)-realizable, so h >= 7 follows from their Theorem 5.2

  Consequently the LOWER bound is theirs. The contribution of this
  paper is the UPPER bound h <= 7 (see h8_impossible.py and the
  Ramsey argument), the uniqueness, and the exact realisation above.

$ python3 h8_impossible.py
======================================================================
The unique (3,4)-Ramsey graph on 8 vertices
======================================================================
  n = 8, m = 16, degrees [4, 4, 4, 4, 4, 4, 4, 4]
  independence number = 2  (need <= 2)   : PASS
  clique number       = 3  (need <= 3)   : PASS
  triangle used as basis: (0, 2, 5)

  12 equations in 9 unknowns over Q
  (non-edges contribute NO equation and NO inequality)

  Groebner basis (grevlex): [1]

  basis is {1}  ->  no solution even over C  : PASS

  Therefore the unique 8-vertex candidate is not realisable,
  and h <= 7.  With h >= 7 (Moser spindle) this gives h = 7.
```
