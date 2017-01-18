TIFFPATH=/nfs/project/PZS0530/setsm/tiff-4.0.3

CC=gcc
CFLAGS=-g -std=c99 -O3 -fopenmp

INCS=-I$(TIFFPATH)/include
LDFLAGS=-L$(TIFFPATH)/lib

setsm : Function.o $(LINKER)
	$(CC) $(CFLAGS) -o setsm Function.o $(LDFLAGS) -lm -ltiff

Function.o : Typedefine.h Function.h Function.c
	$(CC) $(CFLAGS) $(INCS) -c Function.c

.PHONY: clean

clean :
	rm -f setsm
	rm -f *.o

