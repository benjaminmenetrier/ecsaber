ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/vader )

list( APPEND vader_src_files
vader/DefaultCookbook.h
vader/vader.h
vader/recipe.h
vader/recipe.cc
vader/RecipeBase.h
vader/RecipeBase.cc
vader/vader.cc
vader/VaderParameters.h
vader/recipes/AirPotentialTemperature_A.cc
vader/recipes/AirPotentialTemperature_B.cc
vader/recipes/AirPotentialTemperature.h
vader/recipes/AirPressure_A.cc
vader/recipes/AirPressureAtInterface_A.cc
vader/recipes/AirPressureAtInterface_B.cc
vader/recipes/AirPressureAtInterface.h
vader/recipes/AirPressure.h
vader/recipes/AirPressureThickness_A.cc
vader/recipes/AirPressureThickness.h
vader/recipes/AirPressureToKappa_A.cc
vader/recipes/AirPressureToKappa.h
vader/recipes/AirTemperature_A.cc
vader/recipes/AirTemperature_B.cc
vader/recipes/AirTemperature.h
vader/recipes/AirVirtualTemperature_A.cc
vader/recipes/AirVirtualTemperature.h
vader/recipes/CloudIceMixingRatio.h
vader/recipes/CloudIceMixingRatio_A.cc
vader/recipes/CloudLiquidMixingRatio.h
vader/recipes/CloudLiquidMixingRatio_A.cc
vader/recipes/LnAirPressureAtInterface_A.cc
vader/recipes/LnAirPressureAtInterface.h
vader/recipes/DryAirDensityLevelsMinusOne.h
vader/recipes/DryAirDensityLevelsMinusOne_A.cc
vader/recipes/DryAirDensityLevelsMinusOne_B.cc
vader/recipes/HydrostaticExnerLevels.h
vader/recipes/HydrostaticExnerLevels_A.cc
vader/recipes/RainMixingRatio.h
vader/recipes/RainMixingRatio_A.cc
vader/recipes/RelativeHumidity.h
vader/recipes/RelativeHumidity_A.cc
vader/recipes/SurfaceAirPressure.h
vader/recipes/SurfaceAirPressure_A.cc
vader/recipes/TestRecipes.cc
vader/recipes/TestRecipes.h
vader/recipes/TotalMixingRatio.h
vader/recipes/TotalMixingRatio_A.cc
vader/recipes/TotalWater.h
vader/recipes/TotalWater_A.cc
vader/recipes/uwind_at_10m_A.cc
vader/recipes/uwind_at_10m.h
vader/recipes/VirtualPotentialTemperature_A.cc
vader/recipes/VirtualPotentialTemperature_B.cc
vader/recipes/VirtualPotentialTemperature.h
vader/recipes/vwind_at_10m_A.cc
vader/recipes/vwind_at_10m.h
vader/recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater.h
vader/recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater_A.cc
vader/recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater2m.h
vader/recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater2m_A.cc
mo/common_varchange.h
mo/common_varchange.cc
mo/common_linearvarchange.h
mo/common_linearvarchange.cc
mo/constants.h
mo/functions.h
mo/functions.cc
mo/model2geovals_linearvarchange.h
mo/control2analysis_linearvarchange.h
mo/control2analysis_linearvarchange.cc
mo/control2analysis_varchange.h
mo/control2analysis_varchange.cc
mo/eval_air_pressure.h
mo/eval_air_pressure.cc
mo/eval_air_pressure_levels.h
mo/eval_air_pressure_levels.cc
mo/eval_air_temperature.h
mo/eval_air_temperature.cc
mo/eval_cloud_ice_mixing_ratio.h
mo/eval_cloud_ice_mixing_ratio.cc
mo/eval_cloud_liquid_mixing_ratio.h
mo/eval_cloud_liquid_mixing_ratio.cc
mo/eval_dry_air_density.h
mo/eval_dry_air_density.cc
mo/eval_exner.h
mo/eval_exner.cc
mo/eval_geostrophic_to_hydrostatic_pressure.h
mo/eval_geostrophic_to_hydrostatic_pressure.cc
mo/eval_hydrostatic_balance.h
mo/eval_hydrostatic_balance.cc
mo/eval_mio_fields.h
mo/eval_mio_fields.cc
mo/eval_moisture_incrementing_operator.h
mo/eval_moisture_incrementing_operator.cc
mo/eval_moisture_control.h
mo/eval_moisture_control.cc
mo/eval_rain_mixing_ratio.h
mo/eval_rain_mixing_ratio.cc
mo/eval_ratio.h
mo/eval_ratio.cc
mo/eval_relative_humidity.h
mo/eval_relative_humidity.cc
mo/eval_sat_vapour_pressure.h
mo/eval_sat_vapour_pressure.cc
mo/eval_surface_pressure.h
mo/eval_surface_pressure.cc
mo/eval_surface_temperature.h
mo/eval_surface_temperature.cc
mo/eval_surface_wind.h
mo/eval_surface_wind.cc
mo/eval_total_mixing_ratio.h
mo/eval_total_mixing_ratio.cc
mo/eval_total_relative_humidity.h
mo/eval_total_relative_humidity.cc
mo/eval_virtual_potential_temperature.h
mo/eval_virtual_potential_temperature.cc
mo/eval_water_vapor_mixing_ratio.h
mo/eval_water_vapor_mixing_ratio.cc
mo/model2geovals_linearvarchange.h
mo/model2geovals_linearvarchange.cc
mo/model2geovals_varchange.h
mo/model2geovals_varchange.cc
mo/svp_interface.F90
)

if( gsw_FOUND )
  list( APPEND vader_src_files
  vader/recipes/SeaWaterPotentialTemperature.h
  vader/recipes/SeaWaterPotentialTemperature_A.cc
  vader/recipes/SeaWaterTemperature.h
  vader/recipes/SeaWaterTemperature_A.cc
  OceanConversions/OceanConversions.interface.F90
  OceanConversions/OceanConversions.interface.h
  )
endif( gsw_FOUND )

include(GNUInstallDirs)
ecbuild_add_library( TARGET ${PROJECT_NAME}
                     SOURCES ${${PROJECT_NAME}_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
                     LINKER_LANGUAGE CXX )

# Mark Boost as a system library to silence compiler warnings from Boost headers
target_include_directories( ${PROJECT_NAME} SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} )

target_link_libraries( ${PROJECT_NAME} PUBLIC ${oops_LIBRARIES} oops-patch ) #TODO: Change to "oops::oops" once oops adds namespace support

if( gsw_FOUND )
  target_link_libraries( ${PROJECT_NAME} PUBLIC gsw )
  target_compile_definitions( ${PROJECT_NAME} PUBLIC GSW_FOUND )
endif( gsw_FOUND )

#Configure include directory layout for build-tree to match install-tree
set(BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(${PROJECT_NAME}_headers ALL
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_DIR_INCLUDE_PATH}
                    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${BUILD_DIR_INCLUDE_PATH}/${PROJECT_NAME})
add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_headers)

## Include paths
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
