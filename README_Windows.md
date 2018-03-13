# Running SETSM on Windows using Cygwin

Use Cygwin, a large collection of GNU and Open Source tools providing functionality similar to Linux, to build and run SETSM on Windows.

## Installation

1. Download the appropriate installer from https://cygwin.com/install.html
	* setup-x86.exe for 32-bit versions of Windows
	* setup-x86_64.exe for 64-bit versions of Windows
2. Run installer
3. Choose "Install from Internet" as the download source. Then click next.
4. Select root Install directory and then click next.
	* Type in the path to the directory you want Cygwin installed under "Root Directory" (or use the default)
	* Select "All Users (RECOMMENDED)" if Cygwin will be available to all users of the system.
	* Select "Just Me" if you lack Administrator privileges or have specific needs.
5. Type in or browse for the path where you want Cygwin downloaded under "Local Package Directory". Then click next.
6. Select "Direct Connection" for internet connection. Then click next.
7. Choose a download site out of the available sites (does not matter which one you pick). Then click next.
8. Select the following packages by clicking skip and checking Bin (put an X in selected box) for each:
	* Search packages for gcc
		* In Devel category, select to install the gcc-core package (and gcc-g++ package if desired)
	* Search Packages for bin
		* In Devel category, select to install binutils
	* Search packages for make
		* In Devel category, select to install make
		* If wget is required, select to install wget from the Web category
9. Click next to allow for installation

## Install and Run SETSM

1. Download SETSM onto your system in a location you can access it with Cygwin (nearly anywhere on your system).
2. Set up the tiff library on your system as needed following the instructions in the SETSM README.md under Prerequisite.
3. Use the Makefile for gcc (or another compiler if you have downloaded it) to build SETSM according to the SETSM README.md
	* Change the TIFFPATH and uncomment those 3 lines if the tiff library is in a different location.
4. Make sure your images and their corresponding header files (*.xml) have the same filename and are in the same directory.
5. Make sure your default.txt file is in the same directory as where you will run SETSM.
6. Run SETSM according to the SETSM User Manual:

```
./setsm [image1(*.raw|*.tif)] [image2(*.raw|*.tif)] [outputpath/name] [-options]
```

For more SETSM options, refer to the SETSM User Manual or run ./setsm -help

