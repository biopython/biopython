all: ckMeans.so

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $< -o $@

%.so: %module.o
	ld ${LDFLAGS} ${MORELDFLAGS} $< -o $@

%.so: %.o
	ld ${LDFLAGS} ${MORELDFLAGS} $< -o $@

.PHONY: clean
clean: ;
