# Unit chords on the sphere: the extremal almost-equidistant number is 7

Self-contained Overleaf project. Upload this zip, compiler pdfLaTeX, main
document `main.tex`. BibTeX runs automatically.

## The literature check has been done

arXiv:2504.11086 (v2, 22 April 2025) was checked. Result, in short:

* Their paper is titled "Obtuse almost-equiangular sets" and every bound in it
  is for t in [-1, 0]. t = 1/2 is positive and outside all of them. They note
  explicitly that the bound 2(n+1) does NOT extend to positive t, citing the
  Larman-Rogers configuration with alpha(5, 1/5) >= 16 against 2(n+1) = 12.
  So alpha(3, 1/2) is not among their computed values.

* HOWEVER their Theorem 5.2 gives that the 2-Moser spindle MS(2) is
  (3,t)-realizable for t >= t_{2,3} = 0.2553..., and 1/2 exceeds that. MS(2)
  has order 7, with 11 edges and degree sequence [3,3,3,3,3,3,4].
  **The lower bound h >= 7 therefore already follows from their work, and the
  extremal configuration is the Moser spindle.**

The abstract and Section 1 have been rewritten accordingly. The claim is now
the **upper bound h <= 7 and the uniqueness**, which is genuinely outside
everything in the literature, plus the identification of the configuration and
its exact Gram realisation over Q(sqrt 46).

## Where to send it

To the authors of arXiv:2504.11086 — Christine Bachoc (Bordeaux), Bram Bekker
and Fernando de Oliveira Filho (Delft), Philippe Moustrou (Toulouse). Addresses
are on the last page of their paper. They are the right readers: your result
extends their classification to the first positive inner product, and your
Ramsey-plus-Groebner argument is a different technique from their semidefinite
one.

Write to de Oliveira Filho or Bachoc. Keep it short: state that their Theorem
5.2 gives the lower bound, that you have the matching upper bound by a Ramsey
argument plus one Groebner basis, and attach the PDF. Say plainly that you are
an engineer working independently and would be grateful if they checked it.

That framing is the honest one and also the strongest: you are not claiming
their result, you are closing the other side of it.
