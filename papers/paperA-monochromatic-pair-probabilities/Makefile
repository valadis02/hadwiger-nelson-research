MAIN = main
all:
	pdflatex -interaction=nonstopmode $(MAIN).tex
	-bibtex $(MAIN)
	pdflatex -interaction=nonstopmode $(MAIN).tex
	pdflatex -interaction=nonstopmode $(MAIN).tex
clean:
	rm -f *.aux *.log *.out *.bbl *.blg
