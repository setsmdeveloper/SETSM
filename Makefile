TIFFPATH=/fs/project/PZS0530/setsm/tiff-4.0.3

CC=gcc
CFLAGS=-g -std=c99 -O3 -ffast-math -fopenmp -march=native

INCS=-I$(TIFFPATH)/include
LDFLAGS=-L$(TIFFPATH)/lib

setsm : setsm_code.o $(LINKER)
	$(CC) $(CFLAGS) -o setsm setsm_code.o $(LDFLAGS) -lm -ltiff

setsm_code.o : Typedefine.h setsm_code.h setsm_code.c
	$(CC) $(CFLAGS) $(INCS) -c setsm_code.c

.PHONY: clean

clean :
	rm -f setsm
	rm -f *.o

