# Building with cmake on OSC

To build with cmake at OSC:

1. Load the module for your compiler.
2. Ensure that `CXX` is set to that compiler, possibly by running
   `export CXX=$OSC_CXX`.
3. Load the module for libgeotiff: `module load libgeotiff`.
4. Create a directory to do the build in, for example: `cd setsm; mkdir build; cd build`
5. From the build directory, run the following:

```sh
cmake -DCMAKE_PREFIX_PATH="/apps/libgeotiff/1.4.3/tiff;/apps/libgeotiff/1.4.3" /path/to/setsm
```

If your build directory is at `setsm/build`, then this would look like:

```sh
cmake -DCMAKE_PREFIX_PATH="/apps/libgeotiff/1.4.3/tiff;/apps/libgeotiff/1.4.3" ..
```

The `CMAKE_PREFIX_PATH` is a list of paths to search for libraries in. Elements are
separated by semicolons. This tells cmake where to look for the TIFF and GeoTIFF
libraries.
