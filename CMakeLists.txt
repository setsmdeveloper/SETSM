cmake_minimum_required (VERSION 3.0.0)

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

set (GNU_VERSION 4.0.0)
set (INTEL_VERSION 14.0.0)
set (CRAY_VERSION 0.0.0)
set (TIFF_VERSION 4.0.0)
set (GEOTIFF_VERSION 1.4.0)

set (CC_LIST GNU INTEL CRAY)
set (GNU_ID "GNU")
set (INTEL_ID "Intel")
set (CRAY_ID "Cray")

if ($ENV{MPI} MATCHES "on")
    message ("Building with MPI enabled")
    set (GNU_FLAGS -O3 -ffast-math -march=native -DBUILDMPI)
    set (INTEL_FLAGS -O3 -xHost -DBUILDMPI)
    set (CRAY_FLAGS -DBUILDMPI)
else ()
    message ("Building with MPI disabled")
    set (GNU_FLAGS -O3 -ffast-math -march=native)
    set (INTEL_FLAGS -O3 -xHost)
    set (CRAY_FLAGS )
endif()

if (NOT DEFINED OPENMP_OPTION)
    set (OPENMP_OPTION )
endif ()

# Create target
project (setsm CXX)
add_executable (setsm edge_list.cpp grid.cpp grid_triangulation.cpp setsm_code.cpp setsmgeo.cpp)
set_target_properties (setsm PROPERTIES CXX_STANDARD 11)

# Check compiler
foreach (CPL ${CC_LIST})
    if (CMAKE_CXX_COMPILER_ID MATCHES ${${CPL}_ID})
        set (COMPILER ${CPL})
    endif ()
endforeach ()
if (NOT DEFINED COMPILER)
    message (FATAL_ERROR "Unknown compiler ${CMAKE_CXX_COMPILER_ID}.")
endif ()

if (CMAKE_CXX_COMPILER_VERSION LESS ${${COMPILER}_VERSION})
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} minimum version \"${${COMPILER}_VERSION}\" required.")
endif ()

# Find tiff
find_package (TIFF ${TIFF_VERSION})
if (NOT TIFF_FOUND)
    message (FATAL_ERROR "Cannot find tiff/${TIFF_VERSION} or above.")
endif ()

# Find geotiff
find_package (GeoTIFF ${GEOTIFF_VERSION})
if (NOT GEOTIFF_FOUND)
    message (FATAL_ERROR "Cannot find geotiff/${GEOTIFF_VERSION} or above.")
endif ()

# Find OpenMP
find_package (OpenMP)
if (NOT OPENMP_FOUND)
    message (FATAL_ERROR "Cannot find OpenMP.")
endif ()

# Build setsm target
target_compile_options (setsm PUBLIC ${${COMPILER}_FLAGS})
target_include_directories (setsm PUBLIC ${TIFF_INCLUDE_DIR} ${GEOTIFF_INCLUDE_DIR})
target_link_libraries (setsm m ${TIFF_LIBRARIES} ${GEOTIFF_LIBRARIES})

list (APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS}${OPENMP_OPTION})

# Make install
install(TARGETS setsm DESTINATION .)
