###############################################################################
# CMake module to search for PROJ library
#
# On success, the macro sets the following variables:
# PROJ_FOUND       = if the library found
# PROJ_LIBRARY     = full path to the library
# PROJ_INCLUDE_DIR = where to find the library headers 
# also defined, but not for general use are
# PROJ_LIBRARY, where to find the PROJ library.
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
###############################################################################

SET(PROJ_NAMES ${PROJ_NAMES} proj proj_i)

FIND_LIBRARY(PROJ_LIBRARY
    NAMES ${PROJ_NAMES}
    HINTS ENV LD_LIBRARY_PATH) 

get_filename_component(PROJ_LIBRARY_PATH ${PROJ_LIBRARY} DIRECTORY)

FIND_PATH(PROJ_INCLUDE_DIR proj_api.h
  HINTS "${PROJ_LIBRARY_PATH}/../include")

IF(PROJ_INCLUDE_DIR)
  SET(PROJ_VERSION 0)

  SET(PROJ_VERSION_H "${PROJ_INCLUDE_DIR}/proj_api.h")
  FILE(READ ${PROJ_VERSION_H} PROJ_VERSION_H_CONTENTS)

  IF (DEFINED PROJ_VERSION_H_CONTENTS)
    STRING(REGEX REPLACE ".*#define[ \t]PJ_VERSION[ \t]+([0-9]+).*" "\\1" PROJ_VERSION_NUM "${PROJ_VERSION_H_CONTENTS}")

    IF(NOT ${PROJ_VERSION_NULL} MATCHES "[0-9]+")
      MESSAGE(FATAL_ERROR "PROJ version parsing failed!")
    ENDIF()

    IF(PROJ_VERSION_NUM AND NOT "${PROJ_VERSION_NUM}" STREQUAL "0")
      MATH(EXPR PROJ_VERSION_MAJOR "${PROJ_VERSION_NUM} / 100")
      MATH(EXPR PROJ_VERSION_MINOR "${PROJ_VERSION_NUM} % 100 / 10")
      MATH(EXPR PROJ_VERSION_PATCH "${PROJ_VERSION_NUM} % 10")
    ENDIF()

    SET(PROJ_VERSION "${PROJ_VERSION_MAJOR}.${PROJ_VERSION_MINOR}.${PROJ_VERSION_PATCH}"
      CACHE INTERNAL "The version string for PROJ library")

    IF (PROJ_VERSION VERSION_EQUAL PROJ_FIND_VERSION OR
        PROJ_VERSION VERSION_GREATER PROJ_FIND_VERSION)
      MESSAGE(STATUS "Found PROJ version: ${PROJ_VERSION}")
    ELSE()
      MESSAGE(FATAL_ERROR "PROJ version check failed. Version ${PROJ_VERSION} was found, at least version ${PROJ_FIND_VERSION} is required")
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR "Failed to open ${PROJ_VERSION_H} file")
  ENDIF()

ENDIF()



# Handle the QUIETLY and REQUIRED arguments and set PROJ_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PROJ DEFAULT_MSG PROJ_LIBRARY PROJ_INCLUDE_DIR)

IF(PROJ_FOUND)
  SET(PROJ_LIBRARIES ${PROJ_LIBRARY})
ENDIF()
