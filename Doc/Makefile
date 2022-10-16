subdirs :=
BP_VERSION := $(shell python -c "`grep ^__version__ ../Bio/__init__.py`; print('dev' if 'dev' in __version__ else __version__)")

all:  Tutorial.html Tutorial.txt pdf $(subdirs)
pdf:  Tutorial.pdf biopdb_faq.pdf

Tutorial.pdf: Tutorial.tex Tutorial/chapter_*.tex
	pdflatex --shell-escape "\def\bpversion{${BP_VERSION}}\input{Tutorial}"
	pdflatex --shell-escape "\def\bpversion{${BP_VERSION}}\input{Tutorial}"
	pdflatex --shell-escape "\def\bpversion{${BP_VERSION}}\input{Tutorial}"

biopdb_faq.pdf: biopdb_faq.tex
	pdflatex biopdb_faq.tex
	pdflatex biopdb_faq.tex
	pdflatex biopdb_faq.tex

Tutorial.html: Tutorial.tex Tutorial/chapter_*.tex version.sh
	hevea -exec ./version.sh -exec xxdate.exe -fix Tutorial.tex

Tutorial.txt: Tutorial.tex Tutorial/chapter_*.tex version.sh
	hevea -exec ./version.sh -exec xxdate.exe -fix -text Tutorial.tex

version.sh: version.in
	sed "s/\@VERSION\@/${BP_VERSION}/g" version.in > version.sh
	chmod +x version.sh

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
	rm -f Tutorial/*.aux

distclean-subdirs: $(subdirs)
	( for f in $^ ; do $(MAKE) distclean -C $$f ; done )

distclean: clean distclean-subdirs
	rm -f version.sh
	rm -f biopdb_faq.pdf
	rm -f Tutorial.pdf
	rm -f Tutorial.html
	rm -f Tutorial.txt
	rm -f *_motif.gif #output from hacha

.PHONY: $(subdirs)
$(subdirs):
	$(MAKE) -C $@

