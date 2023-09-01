# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Links from SABER
set( targets ErrorCovarianceParameters.h GenericCtlVec.cc GenericCtlVec.h instantiateCovarFactory.h instantiateLocalizationFactory.h )
set( dir src/saber/oops )
foreach( target IN LISTS targets )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_SOURCE_DIR}/ecsaber/saber/${dir}/${target}
                     ${CMAKE_SOURCE_DIR}/ecsaber/saber_patch/${dir}/${target} )
endforeach()

list( APPEND oops_src_files_list
# Parameters
ErrorCovarianceParameters.h

# Drivers
ErrorCovariance.h
Localization.h

# Application
ErrorCovarianceToolbox.h

# Generic control vector
GenericCtlVec.cc
GenericCtlVec.h

# Utilities
Utilities.h
Utilities.cc

# Factories
instantiateCovarFactory.h
instantiateLocalizationFactory.h
)

set( oops_src_files ${oops_src_files_list} PARENT_SCOPE )
