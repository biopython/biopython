# create pdf, txt and html documentation

FILENAME=Installation

all: ${FILENAME}.tex pdf html txt

clean: 
	rm -f ${FILENAME}.aux
	rm -f ${FILENAME}.bbl
	rm -f ${FILENAME}.blg
	rm -f ${FILENAME}.log
	rm -f ${FILENAME}.out
	rm -f ${FILENAME}.toc
	rm -f ${FILENAME}.haux
	rm -f ${FILENAME}.htoc

real_clean: clean
	rm -f ${FILENAME}.pdf
	rm -f ${FILENAME}.txt
	rm -f ${FILENAME}.html

pdf: clean pdf_make
pdf_make: ${FILENAME}.tex
	pdflatex ${FILENAME}.tex
	pdflatex ${FILENAME}.tex
	pdflatex ${FILENAME}.tex

html: clean html_make
html_make: ${FILENAME}.tex
	hevea ${FILENAME}.tex
	hevea ${FILENAME}.tex

txt: clean txt_make
txt_make: ${FILENAME}.tex
	hevea -text ${FILENAME}.tex
	hevea -text ${FILENAME}.tex

