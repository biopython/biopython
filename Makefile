# This is a semi-generic makefile for compiling Python shared
# libraries.
# Created 5/6/2002, jchang

# This should be SOLARIS, IRIX, OSX, or LINUX.
PLATFORM=OSX

# Set this to point to the include files for your current installation
# of python.
PYINCLUDE=/sw/include/python2.2

# Use this to pass any other flags to the C compiler.  This should be
# at least -I$(PYINCLUDE).
BIOPYTHON_ROOT=/Users/jchang/remotecvs/biopython
export BIOPYTHON_ROOT
MORECFLAGS=-I$(PYINCLUDE) -I$(BIOPYTHON_ROOT)/Bio/Tools

# Use this to pass any other flags to the linker.
MORELDFLAGS=

# Specify the files to make.
all: Bio
.PHONY: all Bio manifest sdist

Bio:
	$(MAKE) -C Bio

manifest:
	python setup.py sdist --manifest-only

sdist:
	python setup.py sdist --formats=gztar,zip

# You should not have to change things past here, unless you want to
# tweak flags for your platform.

# Set platform-specific values for these variables.
SOLARIS_CC=gcc
SOLARIS_LD=ld
SOLARIS_CFLAGS=-O2 -Wall
SOLARIS_LDFLAGS=-G

IRIX_CC=gcc
IRIX_LD=ld
IRIX_CFLAGS=-O2 -Wall -g -shared
IRIX_LDFLAGS=-shared -all

OSX_CC=cc
OSX_LD=cc
OSX_CFLAGS=-O2 -Wall -g -traditional-cpp
OSX_LDFLAGS=-bundle -undefined suppress -flat_namespace

LINUX_CC=gcc
LINUX_LD=ld
LINUX_CFLAGS=-O2 -Wall
LINUX_LDFLAGS=-G

CC=$($(PLATFORM)_CC)
LD=$($(PLATFORM)_LD)
CFLAGS=$($(PLATFORM)_CFLAGS)
LDFLAGS=$($(PLATFORM)_LDFLAGS)

# Export these variables to sub-makes.
export PYINCLUDE
export MORECFLAGS
export MORELDFLAGS
export CC
export LD
export CFLAGS
export LDFLAGS

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $< -o $@

%.so: %module.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

%.so: %.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

.PHONY : clean
clean:
	$(MAKE) clean -C Doc
	@find . -name '*.py[oc]' -print -exec rm -f {} \;
	@find . -name '*~' -print -exec rm -f {} \;
	@find . -name '*.so' -print -exec rm -f {} \;
	@find . -name '*.o' -print -exec rm -f {} \;
