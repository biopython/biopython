all: csupport.o cmathfns.so clistfns.so cstringfns.so cdistance.so Classification Clustering
.PHONY: all Classification Clustering

Classification:
	$(MAKE) -C Classification

Clustering:
	$(MAKE) -C Clustering

MORELDFLAGS:=$(MORELDFLAGS) $(BIOPYTHON_ROOT)/Bio/Tools/csupport.o

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $< -o $@

%.so: %module.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

%.so: %.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

.PHONY: clean
clean: ;
