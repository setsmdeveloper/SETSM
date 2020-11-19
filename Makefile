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

COMMON_OBJS = CoordConversion.o SubFunctions.o LSF.o Orthogeneration.o Coregistration.o SDM.o setsmgeo.o grid.o grid_triangulation.o edge_list.o
MPI_OBJS = $(COMMON_OBJS) log_mpi.o
OBJS = $(COMMON_OBJS) log.o
HDRS = Typedefine.hpp CoordConversion.hpp SubFunctions.hpp Template.hpp LSF.hpp Orthogeneration.hpp Coregistration.hpp SDM.hpp setsm_code.hpp setsmgeo.hpp grid_triangulation.hpp grid_types.hpp grid_iterators.hpp basic_topology_types.hpp GridVoxel.hpp git_description.h mpi_helpers.hpp log.hpp

ifeq ($(COMPILER), intel)
  CC=icc
  CXX=icpc
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -O3 -qopenmp -fp-model precise 
  CXXFLAGS=-std=c++11 -O3 -qopenmp -fp-model precise
else ifeq ($(COMPILER), pgi)
  CC=pgcc
  CXX=pgc++
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-c99 -O3 -mp=allcores -fast
  CXXFLAGS=-std=c++11 -O3 -mp=allcores -fast
else ifeq ($(COMPILER), cray)
  CC=cc
  CXX=CC
  MPICC=cc
  MPICXX=CC
  CFLAGS=
  CXXFLAGS=-hstd=c++11 -h aggress
else
  CC=gcc
  CXX=g++
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -g -O3 -fopenmp 
  CXXFLAGS=-std=c++11 -O3 -fopenmp 
endif

$(shell git describe --always --tags --dirty > git_description)
GIT_DESCRIPTION:=$(shell cat git_description)
export GIT_DESCRIPTION

setsm : setsm_code.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o setsm setsm_code.o $(OBJS) $(LDFLAGS) -lm -lgeotiff -ltiff

setsm_mpi : setsm_code_mpi.o $(MPI_OBJS)
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o setsm_mpi setsm_code_mpi.o $(MPI_OBJS) $(LDFLAGS) -lm -lgeotiff -ltiff

setsm_code.o : setsm_code.cpp $(HDRS)
	$(CXX) -c $(CXXFLAGS) $(INCS) setsm_code.cpp -o setsm_code.o

setsm_code_mpi.o : setsm_code.cpp $(HDRS)
	$(MPICXX) -c $(CXXFLAGS) $(MPIFLAGS) $(INCS) setsm_code.cpp -o setsm_code_mpi.o

$(OBJS) : $(HDRS)
$(MPI_OBJS) : $(HDRS)

%.o : %.c
	$(CC) -c $(CFLAGS) $(INCS) $< -o $@

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCS) $< -o $@

%_mpi.o : %.cpp
	$(MPICXX) -c $(CXXFLAGS) $(MPIFLAGS) $(INCS) $< -o $@

.PHONY: clean

clean :
	rm -f setsm setsm_mpi
	rm -f *.o git_description git_description.h

git_description.h: git_description
	echo "#define GIT_DESCRIPTION \"$(GIT_DESCRIPTION)\"" > $@
