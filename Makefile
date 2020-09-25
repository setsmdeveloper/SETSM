# This Makefile requires a compiler-specific Makefile.include.  
# Select the appropriate include file from the Makefile_includes 
# directory, copy to this directory, and rename as Makefile.include.

# If libtiff is installed in a nonstandard location you must edit 
# TIFFPATH and uncomment the following three lines.
TIFFPATH=/projects/sciteam/bazu/setsm/lib/tiff-4.0.3-cray
TIFFINC=-I$(TIFFPATH)/include
TIFFLIB=-L$(TIFFPATH)/lib

# If libgeotiff is installed in a nonstandard location you must edit
# GEOTIFFPATH and uncomment the following three lines.
GEOTIFFPATH=/projects/sciteam/bazu/setsm/lib/geotiff
GEOTIFFINC=-I$(GEOTIFFPATH)/include
GEOTIFFLIB=-L$(GEOTIFFPATH)/lib
PROJLIB=-L/projects/sciteam/bazu/setsm/lib/proj/lib

MPIFLAGS = -DBUILDMPI

INCS = $(TIFFINC) $(GEOTIFFINC)
LDFLAGS = $(TIFFLIB) $(GEOTIFFLIB) $(PROJLIB)

OBJS = CoordConversion.o SubFunctions.o LSF.o Orthogeneration.o Coregistration.o SDM.o setsmgeo.o grid.o grid_triangulation.o edge_list.o log.o
HDRS = Typedefine.hpp CoordConversion.hpp SubFunctions.hpp Template.hpp LSF.hpp Orthogeneration.hpp Coregistration.hpp SDM.hpp setsm_code.hpp setsmgeo.hpp grid_triangulation.hpp grid_types.hpp grid_iterators.hpp basic_topology_types.hpp git_description.h mpi_helpers.hpp log.hpp

ifeq ($(COMPILER), intel)
  CC=icc
  CXX=icpc
  MPICC=mpicc
  MPICXX=mpicxx
  CFLAGS=-std=c99 -O3 -qopenmp -fp-model precise 
  CXXFLAGS=-std=c++11 -O3 -qopenmp -fp-model precise -g
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

$(shell ./update_git_description.sh)
GIT_DESCRIPTION:=$(shell cat git_description)
export GIT_DESCRIPTION

setsm : setsm_code.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o setsm setsm_code.o $(OBJS) $(LDFLAGS) -lm -lgeotiff -ltiff -lz -ljpeg -lproj

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
	rm -f *.o git_description git_description.h

git_description.h: git_description
	echo "#define GIT_DESCRIPTION \"$(GIT_DESCRIPTION)\"" > $@
