PYTHON_INC_DIR=-I/usr/include/python2.3/ -I/usr/include/python2.3/Numeric/
CC=gcc
LEX=flex
 
all:	MMCIFlexmodule.so mmcif_test

MMCIFlexmodule.so: lex.yy.o MMCIFlexmodule.o
	$(CC) -fPIC -shared -o $@ $^ -lfl 

MMCIFlexmodule.o: MMCIFlexmodule.c
	$(CC) -c -fPIC $(PYTHON_INC_DIR) -o $@ $<

#mmcifparser.o: mmcifparser.c
#	$(CC) -c -fPIC $(PYTHON_INC_DIR)  $<

lex.yy.c: mmcif.lex
	$(LEX) mmcif.lex 

lex.yy.o: lex.yy.c
	$(CC) -fPIC -c $<

mmcif_test: mmcif_test.o lex.yy.o
	$(CC) $^ -lfl -o $@

.PHONY:	links clean veryclean

clean:
	rm *.o

veryclean: clean
	- rm *.so lex.yy.c mmcif_test *~
