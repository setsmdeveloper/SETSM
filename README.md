SETSM is a scalable XSEDE system for automated, high resolution terrain 
generation.

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

SETSM is dependent on the TIFF Library, which you must download and install 
separately.  Download LibTIFF from http://libtiff.maptools.org/ and install
it according to the instructions in the README file.  In short,
./configure --prefix=/directory-to-install-in
make
make install

Once LibTIFF is installed you must edit the SETSM Makefile to set the correct
path to the library.  SETSM can then be built simply by typing make.

-------------------------------------------------------------------------------

License

SETSM is released under the Apache 2.0 license, a copy of which is included in
this directory.

