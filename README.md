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

This software includes code derived from the voronoi algorithm by 
Steven Fortune (http://ect.bell-labs.com/who/sjf/) 
as modified by Derek Bradley 
(http://zurich.disneyresearch.com/derekbradley/voronoi.html)

Reference: Steve J. Fortune (1987) A Sweepline Algorithm for Voronoi Diagrams,
Algorithmica 2, 153-174.

-------------------------------------------------------------------------------

Installation instructions

SETSM is dependent on LibTIFF, version 4.0.3 or higher.  Your system may 
already have LibTIFF installed.  If not, you must download and install it 
separately.  Download LibTIFF from http://libtiff.maptools.org/ and install
it according to the instructions in the README file.  In short,
./configure --prefix=/directory-to-install-in
make
make install

Select one of the SETSM Makefile.* files based on the compiler you plan to use.
Copy the selected file to Makefile and edit it if necessary to set the correct 
path to the TIFF library.  SETSM can then be built simply by typing make.

-------------------------------------------------------------------------------

License

SETSM is released under the Apache 2.0 license, a copy of which is included in
this directory.

