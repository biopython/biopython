# create pdf, txt and html documentation

FILENAME=Installation

all: ${FILENAME}.pdf ${FILENAME}.html ${FILENAME}.txt

clean: 
	rm -f ${FILENAME}.aux
	rm -f ${FILENAME}.log
	rm -f ${FILENAME}.out
	rm -f ${FILENAME}.toc
	rm -f ${FILENAME}.haux
	rm -f ${FILENAME}.htoc

distclean: clean
	rm -f ${FILENAME}.pdf
	rm -f ${FILENAME}.txt
	rm -f ${FILENAME}.html

${FILENAME}.pdf: ${FILENAME}.tex
	pdflatex ${FILENAME}.tex
	pdflatex ${FILENAME}.tex
	pdflatex ${FILENAME}.tex

${FILENAME}.html: ${FILENAME}.tex
	hevea -fix ${FILENAME}.tex

${FILENAME}.txt: ${FILENAME}.tex
	hevea -fix -text ${FILENAME}.tex

