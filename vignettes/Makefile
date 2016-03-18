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
	MetabolomicTools.tex\

	
all: 
	make vignette

	
vignette: MetabolomicTools.Rnw
	pdflatex MetabolomicTools.tex
	bibtex MetabolomicTools
	pdflatex MetabolomicTools.tex
	pdflatex MetabolomicTools.tex

	
clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -rf figure/
	rm -rf cache/
