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

### Building SETSM using the provided Makefile

The Makefile supports building SETSM with the Intel, GNU and Cray 
compilers.  Use one of the following commands as appropriate.  If no compiler 
is specified it defaults to gcc/g++.
```
COMPILER=intel make
COMPILER=cray make
make
```

You may need to edit the Makefile to set the 
correct path to the TIFF and GeoTIFF libraries.

You can also set the OPTFLAGS variable to override the default optimization flags,
for example,
```
COMPILER=intel OPTFLAGS='-O3 -fp-model precise' make
```


#### Parallel SETSM with MPI (Message-Passing Interface)
To build SETSM for parallel computing with MPI, follow the above steps then use:
```
make setsm_mpi
```

MPI-parallel SETSM (setsm_mpi) has been tested with MVAPICH2, CrayMPI, and 
IntelMPI.  It should work with other MPI implementations with minor changes 
to the Makefile.

If SETSM is built with both MPI and OpenMP it is usually best to run one 
MPI process per compute node and allow multithreading within the node.

MPI in SETSM **requires** async progress support. Most MPI implementations
provide this, though it may be disabled by default. Usually it can be enabled
using environment variables. To enable async progress for intelmpi, set
`I_MPI_ASYNC_PROGRESS=1`/ Similarly, for mvapich, set `MV2_ASYNC_PROGRESS=1`.

Some MPI implementations do not support async progress. Notably, the MPI
implementation provided by cray has some issues on this front. For these
cases, SETSM can provide it's own async progress support. To enable this,
set `SETSM_CUSTOM_ASYNC=1` in the environment. This mode requires that
the MPI implementation provide multthreaded support. Again, most
implementations support this, but may not have it enabled by default. For
example, some cray systems may require that an environment variable be
set: `MPICH_MAX_THREAD_SAFETY=multiple`. 

#### A note about precision

SETSM can be quite sensitive to roundoff error for some terrains.  We have 
tried to set appropriate compilation flags in the Makefiles, but you may need 
to experiment with floating-point related options if you have problems.

### Building and installing SETSM using CMake

A CMake build is provided as an alternative to the Makefile method described 
above.  The environment variable CXX must be set to the compiler command to be 
used.  If the TIFF library is in a nonstandard location its path must be 
provided using a -D flag as shown in the second example below.
If the TIFF library or the GeoTIFF libraries are in nonstandard locations,
the path must be provided using the `CMAKE_PREFIX_PATH` variable as shown
below.

Example 1:  Build with icpc.

```sh
mkdir build
cd build
CXX=icpc cmake ..
make
make install
```

Example 2:  Build with icpc specifying a location for the TIFF and GeoTIFF
libraries.

```sh
mkdir build
cd build
export CXX=icpc
cmake -DCMAKE_PREFIX_PATH="/apps/libgeotiff/1.4.3/tiff;/apps/libgeotiff/1.4.3" ..
make
make install
```

Example 3: Build with mvapich2 and intel.

```sh
mkdir build
cd build
# below: mpicxx is the intel/mvapich2 compiler
export CXX=mpicxx
cmake \
    -DCMAKE_BUILD_TYPE=mpi \
    -DCMAKE_PREFIX_PATH="/apps/libgeotiff/1.4.3/tiff;/apps/libgeotiff/1.4.3" \
    ..
```


## License

SETSM is released under the Apache 2.0 license, a copy of which is included in
this directory.

