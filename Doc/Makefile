all: doc
.PHONY: all doc clean

doc:
	- pdflatex Tutorial
	- pdflatex Tutorial
	- pdflatex Tutorial
	- pdflatex biopdb_faq.tex 
	- pdflatex biopdb_faq.tex 
	- pdflatex biopdb_faq.tex 
	

clean:
	- rm -f Tutorial.ps Tutorial.log Tutorial.dvi \
	Tutorial.aux 
	- rm -rf api
html:
	- hevea -fix Tutorial.tex
	- hacha Tutorial.html
