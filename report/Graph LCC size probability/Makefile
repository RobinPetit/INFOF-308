report:
	pdflatex report.tex
	bibtex report
	pdflatex report.tex
	pdflatex report.tex
	evince report.pdf &

draft:
	pandoc draft.md -o draft.pdf
	evince draft.pdf &
