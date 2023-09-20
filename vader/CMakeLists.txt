#
# (C) Copyright 2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# VADER_PATCH
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( vader VERSION 1.4.0 LANGUAGES CXX Fortran )

# Links from VADER
links_from_source_repo()

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
include( ${PROJECT_NAME}_compiler_flags )
option( ENABLE_VADER_DOC "Build VADER documentation" OFF )

message( STATUS "VADER variables")
message( STATUS "  - ENABLE_VADER_DOC: ${ENABLE_VADER_DOC}" )

## Dependencies
find_package( Boost 1.64.0 REQUIRED )

# Required
find_package( MPI REQUIRED COMPONENTS C Fortran )
find_package( oops REQUIRED )
find_package( oops-patch REQUIRED )
find_package( NetCDF REQUIRED COMPONENTS C Fortran )

find_package( gsw QUIET )
if( ${gsw_FOUND} )
    message(STATUS "GSW FOUND; Including Marine Recipes gsw (${gsw_VERSION})")
else( ${gsw_FOUND} )
    message(STATUS "GSW NOT FOUND: Excluding Marine Recipes")
endif( ${gsw_FOUND} )


## Sources
add_subdirectory( src )
add_subdirectory( tools )

if( ENABLE_VADER_DOC )
    add_subdirectory( docs )
endif()

## Tests
ecbuild_add_test( TARGET ${PROJECT_NAME}_coding_norms
                  TYPE SCRIPT
                  COMMAND ${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}/src
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()