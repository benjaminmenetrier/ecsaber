# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# ECSABER
################################################################################

# Check oops-jedi status
execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/check_status.sh "oops"
                 RESULT_VARIABLE oops_status )
set_property(GLOBAL PROPERTY oops_status "${oops_status}")

# Check saber-jedi status
execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/check_status.sh "saber"
                 RESULT_VARIABLE saber_status )
set_property(GLOBAL PROPERTY saber_status "${saber_status}")

# Check vader-jedi status
execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/check_status.sh "vader"
                 RESULT_VARIABLE vader_status )
set_property(GLOBAL PROPERTY vader_status "${vader_status}")

function( check_diff projectName srcPath dstPath commandPath )
    get_property(status GLOBAL PROPERTY ${projectName}_status)
#    if( ${status} EQUAL 1 )
        execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff.sh ${srcPath} ${dstPath} ${commandPath} )
#    endif()
endfunction()

function( check_diff_recursively input commandPath )
    get_property(status GLOBAL PROPERTY ${PROJECT_NAME}_status)
#    if( ${status} EQUAL 1 )
      execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/script/check_diff_recursively.sh ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${input} ${commandPath} )
#    endif()
endfunction()
