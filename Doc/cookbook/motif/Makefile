# create pdf, txt and html documentation

FILENAME=motif

all: ${FILENAME}.html ${FILENAME}.pdf ${FILENAME}.txt

clean:
	rm -f ${FILENAME}.aux
	rm -f ${FILENAME}.log
	rm -f ${FILENAME}.out
	rm -f ${FILENAME}.toc
	rm -f ${FILENAME}.haux
	rm -f ${FILENAME}.htoc

distclean: clean
	rm -f ${FILENAME}.pdf
	rm -f ${FILENAME}.html
	rm -f ${FILENAME}.txt

${FILENAME}.pdf: ${FILENAME}.tex
	pdflatex ${FILENAME}.tex
	pdflatex ${FILENAME}.tex
	pdflatex ${FILENAME}.tex

${FILENAME}.html: ${FILENAME}.tex
	hevea -fix ${FILENAME}.tex

${FILENAME}.txt: ${FILENAME}.tex
	hevea -fix -text ${FILENAME}.tex

