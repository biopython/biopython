all: doc
.PHONY: all doc clean

doc:
	- latex Tutorial
	- latex Tutorial
	- latex Tutorial
	- dvips -f Tutorial.dvi > Tutorial.ps

clean:
	rm -f Doc/Tutorial.ps Doc/Tutorial.log Doc/Tutorial.dvi \
	Doc/Tutorial.aux
	rm -rf Doc/api
