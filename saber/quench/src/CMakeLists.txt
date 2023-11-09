# (C) Copyright 2023 Meteorologisk Institutt

list( APPEND quench_src_files

ErrorCovariance.h
Fields.cc
Fields.h
Geometry.cc
Geometry.h
Gom.h
Increment.cc
Increment.h
IncrEnsCtlVec.h
IncrModCtlVec.cc
IncrModCtlVec.h
Interpolator.h
LocalizationMatrix.h
Locations.h
# Logbook.h  # TODO(?): only for latest OOPS versions
Model.h
ModelBias.cc
ModelBias.h
ModelBiasCovariance.h
ModelBiasCtlVec.h
ModelBiasEstimator.h
ModelBiasIncrement.cc
ModelBiasIncrement.h
ObsSpace.h
ObsVec.h
State.cc
State.h
Tlm.h
Variables.h

)

ecbuild_add_library( TARGET   quench
                     SOURCES  ${quench_src_files}
                     PUBLIC_LIBS  oops saber
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quench
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

 target_link_libraries( quench PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )

#Configure include directory layout for build-tree to match install-tree
set(QUENCH_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(quench_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber"
                                     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber/quench)
add_dependencies(quench quench_headers)
target_include_directories(quench PUBLIC $<BUILD_INTERFACE:${QUENCH_BUILD_DIR_INCLUDE_PATH}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>)
