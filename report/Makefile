report:
	latex report
	bibtex report
	latex report
	latex report
	dvips -o report.ps -t letterSize report.dvi
	ps2pdf report.ps report.pdf
	evince report.pdf &

draft:
	pandoc draft.md -o draft.pdf
	evince draft.pdf &

clean:
	rm -f report.aux
	rm -f report.log
	rm -f report.blg
	rm -f report.bbl
	rm -f report.dvi
	rm -f report.ps
	rm -f report.pdf
