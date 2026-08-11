# Hadwiger–Nelson: two computational results

Independent research by **Chrysovalantis Tsirindanis** (National Technical
University of Athens) on the Hadwiger–Nelson problem (the chromatic number of
the plane), carried out within the Polymath16 framework. This repository
collects two papers and the ancillary code that verifies every numerical
claim they make, in exact rational / exact algebraic arithmetic.

## Papers

| | Result | PDF | Source |
|---|---|---|---|
| **Paper A** | First non-trivial quantitative bound for proper 5-colourings of the plane: `p_5(√3) ≥ 1/12`; plus two corrections to the published `k=4` table (`p_4(√7) ≤ 1/3`, `p_4(2) ≥ 2838001/9830400`). | [`papers/paperA-monochromatic-pair-probabilities/paperA.pdf`](papers/paperA-monochromatic-pair-probabilities/paperA.pdf) | [`papers/paperA-monochromatic-pair-probabilities/`](papers/paperA-monochromatic-pair-probabilities/) |
| **Paper B** | The extremal almost-equidistant number on the unit sphere is `h = 7` (upper bound `h ≤ 7` proved via a Gröbner-basis realisability argument; the extremal 7-point configuration is unique and given exactly over `Q(√46)`). | [`papers/paperB-unit-chords-sphere-h7/paperB.pdf`](papers/paperB-unit-chords-sphere-h7/paperB.pdf) | [`papers/paperB-unit-chords-sphere-h7/`](papers/paperB-unit-chords-sphere-h7/) |

Both papers are self-contained LaTeX projects (`main.tex` + `body.tex` +
`shared_preamble.tex` + `refs.bib`); build with `make` (pdflatex → bibtex →
pdflatex ×2) inside each paper's folder.

## Code

`code/` mirrors the "ancillary scripts" appendix of each paper. Every script
is self-contained, re-derives its object from scratch (it does not trust a
cached result), and decides its claim in exact arithmetic; floating point, if
used at all, only proposes candidates.

- **`code/paperA/`** — the linear-programming / exact-hull machinery behind
  Paper A: enumerates proper colourings of the wheel `W_6`
  (`verify_hull.py`), builds the 17- and 19-vertex certificate graphs
  (`build_g17.py`), and independently re-checks the two LP dual certificates
  in `Fraction` arithmetic (`check_certificate.py`,
  `check_certificate19.py`, with their certificate data in
  `certificate.json` / `certificate19.json`; `describe19.py` prints the
  symbolic distance classes of `G_19`).
- **`code/paperB/`** — the sphere computations behind Paper B: the `n ≤ 8`
  Ramsey/Gröbner exclusion (`h_sphere.py`, `h8_impossible.py`), the `n = 7`
  exhaustive search over the 9 isomorphism classes (`h7.py`), and the exact
  Gram-matrix realisation over `Q(√46)` (`h7_gram.py`, `h7_theorem.py`).
- **`code/exploratory/`** — supporting scripts from the same research
  campaign that are not individually cited in either paper's appendix
  (alternative graph-variant builders, an earlier high-precision numerical
  pass, and the original run log). Kept for reproducibility and audit trail.

All scripts were re-run against this repository's copy of the papers before
publication here; every printed assertion matches the corresponding
theorem/table in the PDF.

## Requirements

Python 3.10+, `sympy`. A few scripts also use `numpy`. Install with:

```bash
pip install sympy numpy
```

## Citing

If you use this work, please cite the papers directly (BibTeX entries are not
yet assigned — these are independent, not-yet-peer-reviewed preprints; contact
the author for the current citation form).

## Status

These are **independent, self-published results**, not peer-reviewed
publications. Several bibliography entries in both papers carry an explicit
`CHECK` note where the exact venue, year, or author list of a cited work
could not be confirmed at time of writing — see the References section of
each PDF. Corrections are welcome.
