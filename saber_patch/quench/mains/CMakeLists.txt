# (C) Copyright 2023 Meteorologisk Institutt
 
# VADER libraries
if( vader_FOUND )
    set( vader_LIBRARIES "vader" )
else()
    set( vader_LIBRARIES "" )
endif()

# Executables
ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                        SOURCES quenchConvertState.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_error_covariance_toolbox.x
                        SOURCES quenchErrorCovarianceToolbox.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

