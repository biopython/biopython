all: csupport.o \
	cpairwise2.so \
	cmathfns.so clistfns.so cstringfns.so cdistance.so \
	cSVM.so \
	ckMeans.so

.PHONY: all

MORELDFLAGS:=$(MORELDFLAGS) $(BIOPYTHON_ROOT)/Bio/csupport.o

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $< -o $@

%.so: %module.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

%.so: %.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@



.PHONY: clean
clean: ;
