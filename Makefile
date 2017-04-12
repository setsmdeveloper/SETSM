TIFFPATH=/fs/project/PZS0530/setsm/libs/Owens/tiff-4.0.6

CC=icc
CFLAGS=-g -std=c99 -O3 -qopenmp -xHost

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
