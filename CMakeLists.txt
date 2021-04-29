# minimum version for cmake
cmake_minimum_required (VERSION 3.0.0)

# minimum version for compilers
set (GNU_VERSION 4.9.0) # needed for the openmp 4.5 support
set (INTEL_VERSION 14.0.0)
set (PGI_VERSION 13.0.0)
set (CRAY_VERSION 0.0.0)

# Minimum versions for libraries
set (TIFF_VERSION 4.0.0)
set (GEOTIFF_VERSION 1.4.0)

# The compiler list
set (CC_LIST GNU INTEL PGI CRAY)

# includes geotiff find script
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

# Define compiler IDs
set (GNU_ID "GNU")
set (INTEL_ID "Intel")
set (PGI_ID "PGI")
set (CRAY_ID "Cray")

# MPI build options
# see here: https://stackoverflow.com/a/24470998
#
#TODO it might be nice to do some sanity checking for the user
# if they specify mpi as the CMAKE_BUILD_TYPE, then we might
# want to check if CXX starts with mpi, and if not, print a
# warning
set (CMAKE_CXX_FLAGS_MPI "-DBUILDMPI")

# Set per-compiler flags
set (GNU_FLAGS -O3)
set (INTEL_FLAGS -O3 -fp-model precise)
set (PGI_FLAGS -O3 -mp=allcores -fast)
set (CRAY_FLAGS -hstd=c++11 -h aggress)

# Set standard to c++11 and disable compiler specific extensions
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Create target
project (setsm)

include_directories(${setsm_BINARY_DIR})

# see here: https://stackoverflow.com/a/32062884
add_custom_target(dummy_target ${CMAKE_CURRENT_BINARY_DIR}/git_description.h)

# see header file section here:
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#how-can-i-generate-a-source-file-during-the-build
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/git_description.h
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gen_git_description.sh ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS dummy-command
        )

add_custom_command(
    OUTPUT dummy-command
    COMMAND cmake -E echo
    )

set (SRCS
    CoordConversion.cpp
    SubFunctions.cpp
    LSF.cpp
    Orthogeneration.cpp
    Coregistration.cpp
    SDM.cpp
    setsmgeo.cpp
    grid.cpp
    grid_triangulation.cpp
    edge_list.cpp
    GridVoxel.cpp
    log.cpp
    setsm_code.cpp
)

MESSAGE("CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
add_executable (setsm ${CMAKE_CURRENT_BINARY_DIR}/git_description.h ${SRCS})

# Check compiler
foreach (CPL ${CC_LIST})
	if (CMAKE_CXX_COMPILER_ID MATCHES ${${CPL}_ID})
		set (COMPILER ${CPL})
	endif ()
endforeach (CPL)
if (NOT DEFINED COMPILER)
	message (FATAL_ERROR "Unknown compiler ${CMAKE_CXX_COMPILER_ID}.")
endif (NOT DEFINED COMPILER)
if (CMAKE_CXX_COMPILER_VERSION LESS ${${COMPILER}_VERSION})
	message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} minimum version \"${${COMPILER}_VERSION}\" required.")
endif ()

# Find tiff
find_package (TIFF ${TIFF_VERSION})
if (NOT TIFF_FOUND)
	message (FATAL_ERROR "Cannot find tiff/${TIFF_VERSION} or above.")
endif (NOT TIFF_FOUND)

# Find geotiff
find_package (GeoTIFF)# ${GEOTIFF_VERSION})
if (NOT GEOTIFF_FOUND)
	message (FATAL_ERROR "Cannot find geotiff/${GEOTIFF_VERSION} or above.")
endif (NOT GEOTIFF_FOUND)

# Find OpenMP
find_package (OpenMP)
if (NOT OPENMP_FOUND)
	message (FATAL_ERROR "Cannot find OpenMP.") 
endif (NOT OPENMP_FOUND)

# Build setsm target
target_compile_options (setsm PUBLIC ${${COMPILER}_FLAGS})
target_include_directories (setsm PUBLIC ${TIFF_INCLUDE_DIR})
target_include_directories (setsm PUBLIC ${GEOTIFF_INCLUDE_DIR})
target_link_libraries (setsm m ${TIFF_LIBRARIES})
target_link_libraries (setsm m ${GEOTIFF_LIBRARY})

list (APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})

# Make install
#install(TARGETS setsm DESTINATION .)
