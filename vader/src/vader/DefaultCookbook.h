/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "vader/vader.h"
// Recipe headers
#include "recipes/AirPotentialTemperature.h"
#include "recipes/AirPressureThickness.h"
#include "recipes/AirTemperature.h"
#include "recipes/AirVirtualTemperature.h"
#include "recipes/CloudIceMixingRatio.h"
#include "recipes/CloudLiquidMixingRatio.h"
#include "recipes/DryAirDensityLevelsMinusOne.h"
#include "recipes/HydrostaticExnerLevels.h"
#include "recipes/RainMixingRatio.h"
#include "recipes/TotalMixingRatio.h"
#include "recipes/TotalWater.h"
#include "recipes/uwind_at_10m.h"
#include "recipes/VirtualPotentialTemperature.h"
#include "recipes/vwind_at_10m.h"
#include "recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater.h"
#include "recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater2m.h"

namespace vader {

const cookbookConfigType Vader::defaultCookbookDefinition = {
        // Default VADER cookbook definition
        {"air_temperature",        {AirTemperature_A::Name, AirTemperature_B::Name}},
        {"mass_content_of_cloud_ice_in_atmosphere_layer",
                                   {CloudIceMixingRatio_A::Name}},
        {"mass_content_of_cloud_liquid_water_in_atmosphere_layer",
                                   {CloudLiquidMixingRatio_A::Name}},
        {"dry_air_density_levels_minus_one",  {DryAirDensityLevelsMinusOne_A::Name,
                                               DryAirDensityLevelsMinusOne_B::Name}},
        {"hydrostatic_exner_levels",  {HydrostaticExnerLevels_A::Name}},
        {"potential_temperature",  {AirPotentialTemperature_A::Name}},
        {"qrain",                  {RainMixingRatio_A::Name}},
        {"m_t",                    {TotalMixingRatio_A::Name}},
        {"qt",                     {TotalWater_A::Name}},
        {"specific_humidity",      {WaterVaporMixingRatioWrtMoistAirAndCondensedWater_A::Name}},
        {"specific_humidity_at_two_meters_above_surface",
                                   {WaterVaporMixingRatioWrtMoistAirAndCondensedWater2m_A::Name}},
        {"uwind_at_10m",           {uwind_at_10m_A::Name}},
        {"virtual_potential_temperature",  {VirtualPotentialTemperature_B::Name,
                                            VirtualPotentialTemperature_A::Name}},
        {"virtual_temperature",    {AirVirtualTemperature_A::Name}},
        {"vwind_at_10m",           {vwind_at_10m_A::Name}}
    };

}  // namespace vader
