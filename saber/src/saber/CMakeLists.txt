# (C) Copyright 2017-2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Build list of subdirs with files to add
set( _subdirs blocks bump fastlam generic gsi interpolation oops spectralb util vader )
foreach( _subdir IN LISTS _subdirs )
    add_subdirectory( ${_subdir} )
    list( TRANSFORM ${_subdir}_src_files PREPEND ${_subdir}/ )
    list( APPEND ${PROJECT_NAME}_src_files ${${_subdir}_src_files} )
endforeach()

include( GNUInstallDirs )
ecbuild_add_library( TARGET ${PROJECT_NAME}
                     SOURCES ${${PROJECT_NAME}_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
                     LINKER_LANGUAGE CXX )

if( OpenMP_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX )
endif()

target_link_libraries( ${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )
target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} PUBLIC eckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC fckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC atlas_f )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${oops_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} PUBLIC oops-patch )

if( gsibec_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC gsibec )
    if ( sp_FOUND )
        target_link_libraries( ${PROJECT_NAME} PUBLIC sp::sp_d )
    endif()
endif()
if( FFTW_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC FFTW::fftw3)
endif()
if( vader_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC vader )
endif()
if ( ENABLE_OFFLINE_CODECOV )
    target_link_libraries( ${PROJECT_NAME} PUBLIC gcov )
endif()

# Enable including of local header Fortran file templates
target_include_directories( ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )

# Process fypp files
fckit_target_preprocess_fypp( ${PROJECT_NAME} )

# Fortran module output directory for build and install interfaces
set( MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION} )
set_target_properties( ${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR} )
target_include_directories( ${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}> )
