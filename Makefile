
CC=gcc
OPTIMIZER=-O3
OPENMP=-fopenmp
XHOST= 
LINKER= 


CFLAGS=-g -std=c99 ${OPTIMIZER} ${XHOST} ${OPENMP}

IPATH=-I /nfs/gpfs/PZS0530/setsm/tiff-4.0.3/include
LDPATH=-L /nfs/gpfs/PZS0530/setsm/tiff-4.0.3/lib

setsm : Function.o $(LINKER)
	$(CC) $(CFLAGS) -o setsm Function.o $(LINKER) $(LDPATH) -lm -ltiff

Function.o : Typedefine.h Function.h Function.c
	$(CC) $(CFLAGS) $(IPATH) -c Function.c

stub.o: stub.c
	$(CC) $(CFLAGS) $(IPATH) -c stub.c

.PHONY: clean veryclean

clean :
	rm -f setsm
	rm -f *.o

veryclean: clean
	rm -f core.*
	rm -rf Results/*
	rm -f icc*
	rm -f gcc*
	rm -rf out/	
