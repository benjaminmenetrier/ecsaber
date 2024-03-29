# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# SABER
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( saber VERSION 1.9.0 LANGUAGES CXX Fortran )

# Remove update summary
file( REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )

# Copy and update source from SABER
check_diff_recursively( cmake "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( quench/data "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( src "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( test "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( tools "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( COPYING "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( CPPLINT.cfg "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( LICENSE.md "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( README.md "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
check_diff_recursively( saber-import.cmake.in "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )

# Check update summary
if( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt")
    message( FATAL_ERROR "Merge required:\n bash ${CMAKE_CURRENT_SOURCE_DIR}/update_summary.txt" )
endif()

option(OPENMP "Build saber with OpenMP support" ON)

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
option( ENABLE_SABER_INSTRUMENTATION "Enable SABER instrumentation" OFF )
option( ENABLE_OFFLINE_CODECOV "Enable SABER offline code coverage" OFF )
option( ENABLE_SABER_DOC "Enable SABER documentation" OFF )
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
include( ${PROJECT_NAME}_compiler_flags )

## Dependencies
# Required
if(OPENMP)
  find_package( OpenMP REQUIRED COMPONENTS CXX Fortran )
endif()
find_package( MPI REQUIRED COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS C Fortran )
find_package( eckit 1.24.4 REQUIRED COMPONENTS MPI )
find_package( fckit 0.11.0 REQUIRED )
find_package( atlas 0.35.0 REQUIRED COMPONENTS OMP_Fortran TESSELATION)

if( ENABLE_MKL )
    find_package( MKL )
endif()
if( MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

# OOPS
find_package( oops REQUIRED )
find_package( oops-patch REQUIRED )

# Optionals
find_package( FFTW 3.3.8 QUIET )
find_package( gsibec 1.1.3 QUIET )
find_package( vader 1.6.0 QUIET )

# Optional SABER blocks
if( gsibec_FOUND )
    find_package( sp QUIET )
    message( STATUS "SABER block GSI is enabled" )
else()
    message( STATUS "SABER block GSI is NOT enabled" )
endif()
if( FFTW_FOUND )
    message( STATUS "SABER block FastLAM spectral layer is enabled" )
else()
    message( STATUS "SABER block FastLAM spectral layer is NOT enabled" )
endif()
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    message( STATUS "SABER block SPECTRALB is enabled" )
else()
    message( STATUS "SABER block SPECTRALB is NOT enabled" )
endif()
if( vader_FOUND )
    message( STATUS "SABER blocks based on VADER are enabled" )
else()
    message( STATUS "SABER blocks based on VADER are NOT enabled" )
endif()

## SABER instrumentation
if( ENABLE_SABER_INSTRUMENTATION )
    message( STATUS "SABER instrumentation is enabled" )
    add_definitions(-DENABLE_SABER_INSTRUMENTATION=1)
else()
    message( STATUS "SABER instrumentation is not enabled" )
    add_definitions(-DENABLE_SABER_INSTRUMENTATION=0)
endif()

## Code coverage
if ( ENABLE_OFFLINE_CODECOV )
    message( STATUS "Offline code coverage is enabled" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set( CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -fprofile-arcs -ftest-coverage")
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -fprofile-arcs -ftest-coverage")
else()
    message( STATUS "Offline code coverage is not enabled" )
endif()

## Sources
add_subdirectory( quench )
add_subdirectory( src/saber )
add_subdirectory( tools )
add_subdirectory( test )
if( ENABLE_SABER_DOC )
    message( STATUS "SABER documentation is enabled" )
    add_subdirectory( docs )
else()
    message( STATUS "SABER documentation is not enabled" )
endif()

## Global tests
ecbuild_add_test( TARGET saber_coding_norms_src
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/src )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
