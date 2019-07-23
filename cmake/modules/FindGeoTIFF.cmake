find_path (GEOTIFF_INCLUDE_DIR geotiff.h PATH_SUFFIXES geotiff libgeotiff)

find_library (GEOTIFF_LIBRARIES NAMES geotiff)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GEOTIFF DEFAULT_MSG GEOTIFF_INCLUDE_DIR GEOTIFF_LIBRARIES)
