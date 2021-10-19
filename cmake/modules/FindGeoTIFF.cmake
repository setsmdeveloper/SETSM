###############################################################################
#
# CMake module to search for GeoTIFF library
#
# On success, the macro sets the following variables:
# GEOTIFF_FOUND       = if the library found
# GEOTIFF_LIBRARIES   = full path to the library
# GEOTIFF_INCLUDE_DIR = where to find the library headers
# also defined, but not for general use are
# GEOTIFF_LIBRARY, where to find the PROJ.4 library.
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
# Module source: http://github.com/mloskot/workshop/tree/master/cmake/
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
###############################################################################

SET(GEOTIFF_NAMES geotiff)

IF(WIN32)

    IF(MINGW)
        FIND_PATH(GEOTIFF_INCLUDE_DIR
            geotiff.h
            PATH_PREFIXES geotiff
            PATHS
            /usr/local/include
            /usr/include
            c:/msys/local/include)

        FIND_LIBRARY(GEOTIFF_LIBRARY
            NAMES ${GEOTIFF_NAMES}
            PATHS
            /usr/local/lib
            /usr/lib
            c:/msys/local/lib)
    ENDIF(MINGW)

    IF(MSVC)
        SET(GEOTIFF_INCLUDE_DIR "$ENV{LIB_DIR}/include" CACHE STRING INTERNAL)

        SET(GEOTIFF_NAMES ${GEOTIFF_NAMES} geotiff_i)
        FIND_LIBRARY(GEOTIFF_LIBRARY NAMES
            NAMES ${GEOTIFF_NAMES}
            PATHS
            "$ENV{LIB_DIR}/lib"
            /usr/lib
            c:/msys/local/lib)
    ENDIF(MSVC)

ELSEIF(UNIX)

    FIND_LIBRARY(GEOTIFF_LIBRARY NAMES ${GEOTIFF_NAMES} HINTS ENV LD_LIBRARY_PATH)
    get_filename_component(GEOTIFF_LIBRARY_PATH ${GEOTIFF_LIBRARY} DIRECTORY)

    FIND_PATH(GEOTIFF_INCLUDE_DIR geotiff.h PATH_PREFIXES geotiff HINTS "${GEOTIFF_LIBRARY_PATH}/../include")

ELSE()
    MESSAGE("FindGeoTIFF.cmake: unrecognized or unsupported operating system")
ENDIF()

IF(GEOTIFF_FOUND)
  SET(GEOTIFF_LIBRARIES ${GEOTIFF_LIBRARY})
  message ("SFSDFSDF ${GEOTIFF_LIBRARIES}")
ENDIF()

IF(GEOTIFF_INCLUDE_DIR)
  SET(GEOTIFF_VERSION 0)

  SET(GEOTIFF_VERSION_H "${GEOTIFF_INCLUDE_DIR}/geotiff.h")
  FILE(READ ${GEOTIFF_VERSION_H} GEOTIFF_VERSION_H_CONTENTS)

  IF (DEFINED GEOTIFF_VERSION_H_CONTENTS)
    STRING(REGEX REPLACE ".*#define[ \t]LIBGEOTIFF_VERSION[ \t]+([0-9]+).*" "\\1" GEOTIFF_VERSION_NUM "${GEOTIFF_VERSION_H_CONTENTS}")

    IF(NOT ${GEOTIFF_VERSION_NULL} MATCHES "[0-9]+")
      MESSAGE(FATAL_ERROR "GeoTIFF version parsing failed!")
    ENDIF()

    IF(GEOTIFF_VERSION_NUM AND NOT "${GEOTIFF_VERSION_NUM}" STREQUAL "0")
      MATH(EXPR GTIFF_VERSION_MAJOR "${GEOTIFF_VERSION_NUM} / 1000")
      MATH(EXPR GTIFF_VERSION_MINOR "${GEOTIFF_VERSION_NUM} % 1000 / 100")
      MATH(EXPR GTIFF_VERSION_PATCH "${GEOTIFF_VERSION_NUM} % 100 / 10")
    ENDIF()

    SET(GEOTIFF_VERSION "${GTIFF_VERSION_MAJOR}.${GTIFF_VERSION_MINOR}.${GTIFF_VERSION_PATCH}"
      CACHE INTERNAL "The version string for GeoTIFF library")

    IF (GEOTIFF_VERSION VERSION_EQUAL GeoTIFF_FIND_VERSION OR
        GEOTIFF_VERSION VERSION_GREATER GeoTIFF_FIND_VERSION)
      MESSAGE(STATUS "Found GeoTIFF version: ${GEOTIFF_VERSION}")
    ELSE()
      MESSAGE(FATAL_ERROR "GeoTIFF version check failed. Version ${GEOTIFF_VERSION} was found, at least version ${GeoTIFF_FIND_VERSION} is required")
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR "Failed to open ${GEOTIFF_VERSION_H} file")
  ENDIF()

ENDIF()


# Handle the QUIETLY and REQUIRED arguments and set SPATIALINDEX_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GeoTIFF DEFAULT_MSG GEOTIFF_LIBRARY GEOTIFF_INCLUDE_DIR)
