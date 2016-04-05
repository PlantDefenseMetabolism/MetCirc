LATEXFILES = *.aux\
	*.bbl\
	*.blg\
	*.ilg\
	*.log\
	*.nlo\
	*.nls\
	*.toc\
	*.aux\
	*.dvi\
	*.map\
	*.out\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	MetCirc.tex\

	
all: 
	make vignette

	
vignette: MetCirc.Rnw
	pdflatex MetCirc.tex
	bibtex MetCirc
	pdflatex MetCirc.tex
	pdflatex MetCirc.tex

	
clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -rf figure/
	rm -rf cache/
