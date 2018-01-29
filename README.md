# SETSM

The Surface Extraction from TIN-based Searchspace Minimization (SETSM) software
is a fully automatic algorithm for deriving Digital Terrain Models (DTM) from
pairs of satellite imagery.

The SETSM homepage is [here](https://mjremotesensing.wordpress.com/setsm).

SETSM was written by Myoung-Jong (MJ) Noh, Byrd Polar & Climate Rsch Cntr, 
the [Ohio State University](http://www.osu.edu).
Principal investigator is Ian Howat, School of Earth Sciences,
the [Ohio State University](http://www.osu.edu).
Software engineering team:  Judy Gardiner and Karen Tomko, 
the [Ohio Supercomputer Center](http://www.osc.edu).

This software includes code derived from the voronoi algorithm by 
[Steven Fortune](http://ect.bell-labs.com/who/sjf/) 
as modified by
[Derek Bradley](http://zurich.disneyresearch.com/derekbradley/voronoi.html)

Reference:

 > Steve J. Fortune (1987) A Sweepline Algorithm for Voronoi Diagrams, Algorithmica 2, 153-174.

## Installation instructions

### Prerequisite

SETSM is dependent on LibTIFF, version 4.0.3 or higher.  Your system may 
already have LibTIFF installed.  If not, you must download and install it 
separately.  Download LibTIFF [here](http://libtiff.maptools.org/) and install
it according to the instructions in the README file.  In short,
```
./configure --prefix=/directory-to-install-in
make
make install
```

### Building SETSM using one of the provided Makefiles

Makefiles are provided for building SESTSM with the Intel, PGI, GNU and Cray 
compilers.  Select the appropriate Makefile.* based on the compiler you plan to 
use.  Copy the selected file to Makefile and edit it if necessary to set the 
correct path to the TIFF library.  SETSM can then be built simply by typing:
```
make
```

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
```

Example 2:  Build with icc specifying a location for the TIFF library.
```
mkdir build
cd build
CC=icc cmake -DCMAKE_PREFIX_PATH=$HOME/libTIFF ..
make
```

Note:  The SETSM CMake build with the Cray compiler does not yet work.  

## License

SETSM is released under the Apache 2.0 license, a copy of which is included in
this directory.

