# (C) Copyright 2023 Meteorologisk Institutt

list( APPEND oops_modified_src_files
base/ModelSpaceCovarianceParametersBase.h
interface/ModelData.h
mpi/mpi.cc
mpi/mpi.h
runs/ConvertState.h
)

list( APPEND oops_jedi_src_files
base/DataSetBase.h
base/FieldSet3D.cc
base/FieldSet3D.h
base/FieldSet4D.cc
base/FieldSet4D.h
base/ParameterTraitsVariables.cc
base/ParameterTraitsVariables.h
base/WriteParametersBase.h
util/AnyOf.h
util/AssociativeContainers.h
util/AtlasArrayUtil.cc
util/AtlasArrayUtil.h
util/CompositePath.cc
util/CompositePath.h
util/ConfigFunctions.cc
util/ConfigFunctions.h
util/FieldSetHelpers.cc
util/FieldSetHelpers.h
util/FieldSetOperations.cc
util/FieldSetOperations.h
util/FloatCompare.h
util/gatherPrint.h
util/NamedEnumerator.h
util/parameters/ConfigurationParameter.cc
util/parameters/ConfigurationParameter.h
util/parameters/EmptyParameters.h
util/parameters/GenericParameters.h
util/parameters/HasDiracParameters_.h
util/parameters/HasParameters_.h
util/parameters/HasReadParameters_.h
util/parameters/HasWriteParameters_.h
util/parameters/IgnoreOtherParameters.cc
util/parameters/IgnoreOtherParameters.h
util/parameters/NumericConstraints.h
util/parameters/ObjectJsonSchema.cc
util/parameters/ObjectJsonSchema.h
util/parameters/OptionalParameter.h
util/parameters/OptionalParameter.cc
util/parameters/OptionalPolymorphicParameter.h
util/parameters/Parameter.h
util/parameters/ParameterBase.cc
util/parameters/ParameterBase.h
util/parameters/ParameterConstraint.h
util/parameters/Parameters.cc
util/parameters/Parameters.h
util/parameters/ParametersOrConfiguration.h
util/parameters/ParameterTraits.cc
util/parameters/ParameterTraits.h
util/parameters/ParameterTraitsAnyOf.h
util/parameters/ParameterTraitsScalarOrMap.h
util/parameters/PolymorphicParameter.h
util/parameters/PolymorphicParameterTraits.h
util/parameters/PropertyJsonSchema.cc
util/parameters/PropertyJsonSchema.h
util/parameters/RequiredParameter.h
util/parameters/RequiredPolymorphicParameter.h
util/PartialDateTime.cc
util/PartialDateTime.h
util/Random.h
util/ScalarOrMap.h
util/stringFunctions.cc
util/stringFunctions.h
util/TypeTraits.h
)

file( MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/util )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/util/parameters )

foreach( src_file ${oops_jedi_src_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_CURRENT_SOURCE_DIR}/../../../oops-jedi/src/oops/${src_file}
                     ${CMAKE_CURRENT_SOURCE_DIR}/${src_file} )
endforeach()

list( APPEND oops_ecmwf_src_files
util/abor1_cpp.h
util/dateFunctions.cc
util/dateFunctions.h
util/DateTime.cc
util/DateTime.h
util/Duration.cc
util/Duration.h
util/formats.h
util/IntSetParser.cc
util/IntSetParser.h
util/kinds.F90
util/Logger.h
util/missingValues.cc
util/missingValues.h
util/ObjectCounter.h
util/Printable.h
util/Timer.h
)

foreach( src_file ${oops_ecmwf_src_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_SOURCE_DIR}/oops/src/${src_file}
                     ${CMAKE_CURRENT_SOURCE_DIR}/${src_file} )
endforeach()

list( APPEND oops-patch_src_files
${oops_modified_src_files}
${oops_jedi_src_files}
${oops_ecmwf_src_files}
)

include(GNUInstallDirs)
ecbuild_add_library( TARGET oops-patch
                     SOURCES ${oops-patch_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oops-patch
                     LINKER_LANGUAGE CXX )

target_compile_features( oops-patch PUBLIC cxx_std_14 )
if(OpenMP_FOUND)
    target_link_libraries( oops-patch PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX )
endif()

target_link_libraries( oops-patch PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )
target_link_libraries( oops-patch PUBLIC Eigen3::Eigen )
target_link_libraries( oops-patch PUBLIC eckit )
target_link_libraries( oops-patch PUBLIC fckit )
target_link_libraries( oops-patch PUBLIC atlas_f )
target_link_libraries( oops-patch PUBLIC Boost::boost )
target_link_libraries( oops-patch PUBLIC oops )
target_link_libraries( oops-patch PUBLIC util-patch )

# Enable including of local header Fortran file templates
target_include_directories( oops-patch PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )

# Fortran module output directory for build and install interfaces
set( MODULE_DIR module/oops-patch/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION} )
set_target_properties( oops-patch PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR} )
target_include_directories( oops-patch INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}> )
