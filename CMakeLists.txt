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

# Check oops-ecmwf additions from oops-jedi
execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/check_oops.sh )

# Check saber-jedi status
execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/check_status.sh "saber"
                 RESULT_VARIABLE saber_status )
set_property(GLOBAL PROPERTY saber_status "${saber_status}")

# Check vader-jedi status
execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/check_status.sh "vader"
                 RESULT_VARIABLE vader_status )
set_property(GLOBAL PROPERTY vader_status "${vader_status}")

function( copy_and_update_recursively input )
    get_property(status GLOBAL PROPERTY ${PROJECT_NAME}_status)
    if( ${status} EQUAL 1 )
      execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/copy_and_update_recursively.sh ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${input} )
    endif()
endfunction()

function( copy_and_update input output)
    get_property(oops_status GLOBAL PROPERTY oops_status)
    if( ${oops_status} EQUAL 1 )
      execute_process( COMMAND ${CMAKE_SOURCE_DIR}/ecsaber/copy_and_update.sh ${input} ${output})
    endif()
endfunction()
