all: cfastpairwise.so csupport.so
.PHONY: all

cfastpairwise.so: cfastpairwisemodule.o
	$(LD) $(LDFLAGS) $(MORELDFLAGS) $(BIOPYTHON_ROOT)/Bio/Tools/csupport.o $< -o $@

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $< -o $@

%.so: %module.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

%.so: %.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

.PHONY: clean
clean: ;
