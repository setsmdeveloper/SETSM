# This Makefile requires a compiler-specific Makefile.include.  
# Select the appropriate include file from the Makefile_includes 
# directory, copy to this directory, and rename as Makefile.include.

# If libtiff is installed in a nonstandard location you must edit 
# TIFFPATH and uncomment the following three lines.
#TIFFPATH=$(HOME)/libtiff
#TIFFINC=-I$(TIFFPATH)/include
#TIFFLIB=-L$(TIFFPATH)/lib

# If libgeotiff is installed in a nonstandard location you must edit
# GEOTIFFPATH and uncomment the following three lines.
#GEOTIFFPATH=$(SETSMHOME)/libgeotiff-1.4.2
#GEOTIFFINC=-I$(GEOTIFFPATH)/include
#GEOTIFFLIB=-L$(GEOTIFFPATH)/lib

MPIFLAGS = -DBUILDMPI

INCS = $(TIFFINC) $(GEOTIFFINC)
LDFLAGS = $(TIFFLIB) $(GEOTIFFLIB)

OBJS = setsmgeo.o voronoi_setsm.o
HDRS = Typedefine.hpp setsm_code.hpp setsmgeo.hpp voronoi_setsm.hpp


ifeq ($(COMPILER), intel)
  CC=icc
  CXX=icpc
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -O3 -qopenmp -xHost -fp-model precise
  CXXFLAGS=-O3 -qopenmp -xHost -fp-model precise
else ifeq ($(COMPILER), pgi)
  CC=pgcc
  CXX=pgc++
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-c99 -O3 -mp=allcores -fast
  CXXFLAGS=-O3 -mp=allcores -fast
else
  CC=gcc
  CXX=g++
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -O3 -fopenmp -march=native -ffast-math
  CXXFLAGS=-O3 -fopenmp -march=native -ffast-math
endif

setsm : setsm_code.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o setsm setsm_code.o $(OBJS) $(LDFLAGS) -lm -lgeotiff -ltiff

setsm_mpi : setsm_code_mpi.o $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o setsm_mpi setsm_code_mpi.o $(OBJS) $(LDFLAGS) -lm -lgeotiff -ltiff

setsm_code.o : setsm_code.cpp $(HDRS)
	$(CXX) -c $(CXXFLAGS) $(INCS) setsm_code.cpp -o setsm_code.o

setsm_code_mpi.o : setsm_code.cpp $(HDRS)
	$(MPICXX) -c $(CXXFLAGS) $(MPIFLAGS) $(INCS) setsm_code.cpp -o setsm_code_mpi.o

$(OBJS) : $(HDRS)

%.o : %.c
	$(CC) -c $(CFLAGS) $(INCS) $< -o $@

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCS) $< -o $@

.PHONY: clean

clean :
	rm -f setsm setsm_mpi
	rm -f *.o

