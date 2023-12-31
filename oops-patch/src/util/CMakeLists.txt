# (C) Copyright 2023 Meteorologisk Institutt

list( APPEND util-patch_src_files
LibOOPS.h
)

foreach( src_file ${util-patch_src_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_SOURCE_DIR}/oops/src/util/${src_file}
                     ${CMAKE_CURRENT_SOURCE_DIR}/${src_file} )
endforeach()

include(GNUInstallDirs)
ecbuild_add_library( TARGET util-patch
                     SOURCES ${util-patch_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/util-patch
                     LINKER_LANGUAGE CXX )

if(OpenMP_FOUND)
    target_link_libraries( util-patch PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX )
endif()
target_link_libraries( util-patch PUBLIC eckit )

# Enable including of local header Fortran file templates
target_include_directories( util-patch PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )
