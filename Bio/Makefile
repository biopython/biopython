all: csupport.o \
	cpairwise2.so \
	cmathfns.so clistfns.so cstringfns.so cdistance.so \
	cMarkovModel.so \
	trie.so \
	cSVM.so \
	ckMeans.so

.PHONY: all

trie: trie.c
	${CC} -o trie trie.c

trie.so: triemodule.o trie.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} trie.o $< -o $@

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
