# Hadwiger–Nelson: results and verification code

Independent research by **Chrysovalantis Tsirindanis** (National Technical University of Athens) on the
Hadwiger–Nelson problem (the chromatic number of the plane) and related questions, carried out in the
Polymath16 framework. The repository holds every paper and the code that re-checks each computational
claim in them. All checks use exact arithmetic (rationals, or exact algebraic numbers). Floating point,
where it appears at all, only proposes candidates.

These are preprints and have not been peer reviewed.

## Results

| | Result | Paper |
|---|---|---|
| A | `p_5(√3) ≥ 1/12` for proper 5-colourings of the plane (first non-trivial bound); corrections `p_4(√7) ≤ 1/3`, `p_4(2) ≥ 2838001/9830400` to the published `k = 4` table | [paperA.pdf](papers/paperA-monochromatic-pair-probabilities/paperA.pdf) |
| B | Unit chords on the sphere: extremal almost-equidistant number `h = 7`, with the extremal set unique and exact over `Q(√46)` | [paperB.pdf](papers/paperB-unit-chords-sphere-h7/paperB.pdf) |
| C | Forbidden distances `{1, d}`, `d = (√2+√6)/2`: **`7 ≤ χ({1,d}) ≤ 12`**. Also an explicit vertex-critical 6-chromatic 30-vertex `{1,d}`-graph `G30`, a mod-8 obstruction proposition, and the limits of the method (no `≥ 8`, no `χ(ℝ²) ≥ 6` this way) | [paperC.pdf](papers/paperC-two-distance-chi7/paperC.pdf) |
| D | Map-type colourings: Sokolov–Voronov obstruction radius is exactly `ρ + 2`, `ρ ≤ 1/√3`; a parity proof of their "eleventh trichromatic point" step, plus an uncovered case recorded as open | [paperD.pdf](papers/paperD-map-type-obstruction-radius/paperD.pdf) |
| all | All of the above in one document | [combined.pdf](papers/combined/combined.pdf) |

### Attribution for result C

The `≥ 7` proof builds on Polymath16 work: the gadget `K8−4` is due to Parts (thread 15, 24/12/2019);
the observation that `p = ½` forces a monochromatic leg of a rigid triangle is de Grey's (28/12/2019);
the strict averaging lemma `p_δ ≤ ½` was proved by Ágoston for `k = 4` and sketched for all `k` by
Pálvölgyi (thread 12). The finite SAT certificate, `G30`, the upper bound and the obstruction and limit
results are new here. Earlier `{1,d}` lower bounds: φ (Huddleston, via Parts' 16-vertex graph), `d = 2`
(Exoo–Ismailescu).

## Verifying

```bash
pip install sympy numpy mpmath python-sat     # plus gcc
./verify_all.sh            # ~3 minutes
./verify_all.sh --full     # also the slow Paper B searches (~9 minutes more)
```

GitHub Actions runs `verify_all.sh` on every push (`.github/workflows/verify.yml`; start a run by hand
with `full` ticked for the long checks).

| Check | Script | What it decides |
|---|---|---|
| A | `code/paperA/*.py` | `W6` hull, the `G17`/`G19` graphs, both LP dual certificates in `Fraction` arithmetic |
| B | `code/paperB/h8_impossible.py`, `h7_theorem.py` (`--full`: `h7.py`, `h_sphere.py`, `h7_gram.py`) | `n = 8` not realisable; exact 7-point configuration |
| C, ≤ 12 | `code/two_distance/verify_upper12.py` | 12-colouring sublattice has index 12 and separation exactly 2 × tile diameter |
| C, G30 | `code/two_distance/verify_g30.py` | edges are exactly the pairs at distance 1 or `d`; 6-colouring proper; no 5-colouring (DRUP proof, checked); 30/30 vertex deletions 5-colourable |
| C, ≥ 7 | `code/two_distance/prove_chi7.sh` | rebuilds the 262-point configuration from `cert7.json` in `Q(i,√3,√7)`; checks the gadget, the 1734 edges and the 248 rigid triangles; CaDiCaL refutes 6-colourability; the ~600k-step DRUP proof is checked by the bundled `rup_checker.c` (and by `drat-trim` if it is on `PATH`) |

`code/two_distance/qf.py` is the exact arithmetic used by the C checks: elements of `Q(√3,√7)` as
4-tuples of `Fraction`s over the basis `1, √3, √7, √21`.

What a machine check covers for C: the finite combinatorial core. The analytic part (amenability of
`E(2)`, the invariant mean, the averaging lemma) is in the paper, not in code.

## Layout

```
papers/   paperA-*, paperB-* (LaTeX, build with make), paperC-two-distance-chi7 (LaTeX),
          paperD-map-type-obstruction-radius (PDF only), combined (LaTeX)
code/     paperA/, paperB/, two_distance/   verification scripts
          exploratory/, two_distance/exploratory/   research scripts kept as an audit trail
```

## Known gaps

- Paper D: only the PDF is here. The LaTeX source and its verification scripts are not yet in the repository.
- `code/two_distance/exploratory/` was used to find the results, not to check them. Some scripts need
  PARI/GP or files they produced at the time. They are not run by `verify_all.sh`.
- `code/paperA/describe19.py` imports `bigbuild2` from `code/exploratory/` (`verify_all.sh --full`
  sets `PYTHONPATH`).
- No Lean formalisation yet. A partial formalisation
  of the finite lemmas is possible future work.
