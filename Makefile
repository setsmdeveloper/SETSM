# This Makefile requires a compiler-specific Makefile.include.  
# Select the appropriate include file from the Makefile_includes 
# directory, copy to this directory, and rename as Makefile.include.

# If libtiff is installed in a nonstandard location you must edit 
# TIFFPATH and uncomment the following three lines.
TIFFPATH=/home/noh.56/software/tiff-4.0.3/libtiff
TIFFINC=-I/home/noh.56/software/tiff-4.0.3/include
TIFFLIB=-L/home/noh.56/software/tiff-4.0.3/lib

# If libgeotiff is installed in a nonstandard location you must edit
# GEOTIFFPATH and uncomment the following three lines.
GEOTIFFPATH=/home/noh.56/software/libgeotiff-1.4.2/libxtiff
GEOTIFFINC=-I/home/noh.56/software/libgeotiff-1.4.2/include
GEOTIFFLIB=-L/home/noh.56/software/libgeotiff-1.4.2/lib

MPIFLAGS = -DBUILDMPI

INCS = $(TIFFINC) $(GEOTIFFINC)
LDFLAGS = $(TIFFLIB) $(GEOTIFFLIB)

OBJS = setsmgeo.o grid.o grid_triangulation.o edge_list.o
HDRS = Typedefine.hpp setsm_code.hpp setsmgeo.hpp grid_triangulation.hpp grid_types.hpp grid_iterators.hpp basic_topology_types.hpp


ifeq ($(COMPILER), intel)
  CC=icc
  CXX=icpc
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -qopenmp -fp-model precise 
  CXXFLAGS=-g -std=c++11 -qopenmp -fp-model precise
else ifeq ($(COMPILER), pgi)
  CC=pgcc
  CXX=pgc++
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-c99 -O3 -mp=allcores -fast
  CXXFLAGS=-std=c++11 -O3 -mp=allcores -fast
else
  CC=gcc
  CXX=g++
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -O3 -fopenmp 
  CXXFLAGS=-O3 -fopenmp 
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

