# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# ECSABER
################################################################################

function( check_diff projectName srcPath dstPath commandPath )
    execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff.sh ${srcPath} ${dstPath} ${commandPath} )
endfunction()

function( check_diff_recursively input commandPath )
    execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff_recursively.sh ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${input} ${commandPath} )
endfunction()
