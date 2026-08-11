# Monochromatic-pair probabilities for 5-colourings of the plane

Self-contained Overleaf project. Upload this zip, compiler pdfLaTeX, main
document `main.tex`. BibTeX runs automatically.

## Where to send it

**Polymath16 research thread 17**, on Dustin Mixon's blog. It is still open,
Mixon has explicitly invited comments, and the paper answers a question de Grey
posed there directly: he wrote that for k >= 5 the only facts on record are
p_1 = 0 and p_d < 1/2 for d > 1/2. That is the gap this fills.

A public comment beats a cold email: it is read by everyone involved at once,
and if something is wrong you will be told politely rather than ignored.

Keep the comment to two paragraphs — what you found (the 1/12 bound and the two
corrections to the k=4 table), a link to the PDF, a link to the scripts, and one
line saying you are independent and would welcome checking.

## Before sending

Put the ancillary scripts in a public repository and link them in the paper.
`check_certificate19.py` is the strongest single item in the submission: it
rebuilds the graph, re-derives the profiles, and verifies the rational
inequalities without trusting any solver. A reader who can re-run it in one
command will trust the rest.

## The question you will be asked first

"Is 0.2886964 actually new, or did someone already have it?" That is why the
Attribution section is at the front rather than in a footnote. Have the answer
ready and the conversation starts from the mathematics.
