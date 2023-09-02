# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# ECSABER
################################################################################

function(links_from_source_repo)
    file( STRINGS .gitignore targets )
    string( REPLACE "${PROJECT_NAME}_patch" "${PROJECT_NAME}" REPO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
    foreach( target IN LISTS targets )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${REPO_SOURCE_DIR}/${target}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${target} )
    endforeach()
endfunction()

add_subdirectory( oops_patch )
add_subdirectory( vader_patch )
add_subdirectory( saber_patch )
