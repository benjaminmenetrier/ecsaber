# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Links from SABER
set( targets AtlasInterpWrapper.cc AtlasInterpWrapper.h VertProj.cc VertProj.h GaussToCS.cc GaussToCS.h )
set( dir src/saber/interpolation )
foreach( target IN LISTS targets )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_SOURCE_DIR}/ecsaber/saber/${dir}/${target}
                     ${CMAKE_SOURCE_DIR}/ecsaber/saber_patch/${dir}/${target} )
endforeach()

list(APPEND interpolation_src_files_list
# ATLAS interpolation wrapper
AtlasInterpWrapper.cc
AtlasInterpWrapper.h

# Interpolation block
#Geometry.cc
#Geometry.h
#Interpolation.cc
#Interpolation.h

# VertProj block
VertProj.cc
VertProj.h
)

if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    list(APPEND interpolation_src_files_list
        # ATLAS gauss to cubed-sphere-dual saber block
        GaussToCS.cc
        GaussToCS.h
    )
endif()

set( interpolation_src_files

${interpolation_src_files_list}

PARENT_SCOPE
)
