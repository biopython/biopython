all: cSVM.so

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $< -o $@

%.so: %module.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

%.so: %.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $< -o $@

.PHONY: clean
clean: ;
