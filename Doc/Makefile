subdirs := install


all:  Tutorial.html Tutorial.txt pdf $(subdirs)
pdf:  Tutorial.pdf biopdb_faq.pdf

Tutorial.pdf: Tutorial.tex
	pdflatex Tutorial.tex
	pdflatex Tutorial.tex
	pdflatex Tutorial.tex

biopdb_faq.pdf: biopdb_faq.tex
	pdflatex biopdb_faq.tex
	pdflatex biopdb_faq.tex
	pdflatex biopdb_faq.tex

Tutorial.html: Tutorial.tex
	hevea -fix Tutorial.tex

Tutorial.txt: Tutorial.tex
	hevea -fix -text Tutorial.tex

clean-subdirs: $(subdirs)
	( for f in $^ ; do $(MAKE) clean -C $$f ; done )

clean: clean-subdirs
	rm -f Tutorial.aux
	rm -f Tutorial.toc
	rm -f Tutorial.log
	rm -f Tutorial.out
	rm -f Tutorial.haux
	rm -f Tutorial.htoc
	rm -f biopdb_faq.aux
	rm -f biopdb_faq.log
	rm -f biopdb_faq.out

distclean-subdirs: $(subdirs)
	( for f in $^ ; do $(MAKE) distclean -C $$f ; done )

distclean: clean distclean-subdirs
	rm -f biopdb_faq.pdf
	rm -f Tutorial.pdf
	rm -f Tutorial.html
	rm -f Tutorial.txt
	rm -f *_motif.gif #output from hacha

.PHONY: $(subdirs)
$(subdirs):
	$(MAKE) -C $@

