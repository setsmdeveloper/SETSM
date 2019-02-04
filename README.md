# SETSM

The Surface Extraction from TIN-based Searchspace Minimization (SETSM) software
is a fully automatic algorithm for deriving Digital Terrain Models (DTM) from
pairs of satellite imagery.

SETSM homepage, mjremotesensing.wordpress.com/setsm

SETSM was written by Myoung-Jong (MJ) Noh, Byrd Polar & Climate Rsch Cntr, 
the Ohio State University, www.osu.edu.
Principal investigator is Ian Howat, School of Earth Sciences,
the Ohio State University, www.osu.edu.
Software engineering team:  Judy Gardiner and Karen Tomko, 
the Ohio Supercomputer Center, www.osc.edu

The triangulation algorithm in the code was adapted from the methods in these 
two papers:
[1] Wenzhou Wu, Yikang Rui, Fenzhen Su, Liang Cheng & Jiechen Wang (2014) 
Novel parallel algorithm for constructing Delaunay triangulation based on a
twofold-divide-and-conquer scheme, GIScience & Remote Sensing, 51:5, 537-554, DOI:
10.1080/15481603.2014.946666
[2] Leonidas Guibas & Jorge Stolfi (1985) Primitives for the Manipulation of 
General Subdivisions and the Computation of Voronoi Diagrams, ACM Transactions on 
Graphics, 51:2, 74-123.

## Installation instructions

### Prerequisites
SETSM is dependent on LibTIFF, version 4.0.3 or higher.  Your system may 
already have LibTIFF installed.  If not, you must download and install it 
separately.  Download LibTIFF from http://libtiff.maptools.org/ and install
it according to the instructions in the README file.  In short,
```
./configure --prefix=/directory-to-install-in
make
make install
```
SETSM is also dependent on LibGeoTIFF, version 1.4.0 or higher.  Your system
may already have LibGeoTIFF isntalled. If not, you must download and install
LibTIFF (see above) and PROJ, as they are dependencies of LibGeoTIFF.  PROJ
can be downloaded from https://download.osgeo.org/proj/ and installed according
to the instructions in the README file.  After both LibTIFF and PROJ are installed,
LibGeoTIFF can be downloaded from https://download.osgeo.org/geotiff/libgeotiff/
and installed according to the instructions in the README file. In short,
```
./configure --prefix=/directory-to-install-in \
--with-proj=/directory-of-proj-build \
--with-libtiff=/directory-of-libtiff-build
make
make install
```

### Building SETSM using one of the provided Makefiles

Makefiles are provided for building SETSM with the Intel, PGI, GNU and Cray 
compilers.  Select the appropriate Makefile.* based on the compiler you plan to 
use.  Copy the selected file to Makefile and edit it if necessary to set the 
correct path to the TIFF and GeoTIFF libraries.  SETSM can then be built simply by typing:
```
make
```

#### Parallel SETSM with MPI (Message-Passing Interface)
To build SETSM for parallel computing with MPI, follow the above steps then use:
```
make setsm_mpi
```

MPI-parallel SETSM (setsm_mpi) has been tested with MVAPICH2 and 
OpenMPI.  It should work with other MPI implementations with minor changes 
to the Makefile.

If SETSM is built with both MPI and OpenMP it is usually best to run one 
MPI process per compute node and allow threading within the node.

#### A note about precision

SETSM can be quite sensitive to roundoff error for some terrains.  We have 
tried to set appropriate compilation flags in the Makefiles, but you may need 
to experiment with floating-point related options if you have problems.

### Building and installing SETSM using CMake

A CMake build is provided as an alternative to the Makefile method described 
above.  The environment variable CC must be set to the compiler command to be 
used.  If the TIFF library is in a nonstandard location its path must be 
provided using a -D flag as shown in the second example below.

Example 1:  Build with icc.
```
mkdir build
cd build
CC=icc cmake ..
make
make install
```

Example 2:  Build with icc specifying a location for the TIFF library.
```
mkdir build
cd build
CC=icc cmake -DCMAKE_PREFIX_PATH=$HOME/libTIFF ..
make
make install
```

Example 3:  Build with icc specifying a location for make install.
```
mkdir build
cd build
CC=icc cmake -DCMAKE_INSTALL_PREFIX=./desired_location ..
make
make install
```

Example 4: Build with mvapich2.
```
mkdir build
cd build
CC=mpicc MPI=on cmake ..
make
make install
```

Note:  The SETSM CMake build with the Cray compiler does not yet work.  

## License

SETSM is released under the Apache 2.0 license, a copy of which is included in
this directory.

