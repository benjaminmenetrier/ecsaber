/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "oops/util/Logger.h"
#include "vader/recipes/SeaWaterPotentialTemperature.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char SeaWaterPotentialTemperature_A::Name[] = "SeaWaterPotentialTemperature_A";
const std::vector<std::string> SeaWaterPotentialTemperature_A::Ingredients = {
                                                                "sea_water_temperature",
                                                                "sea_water_salinity",
                                                                "latitude",
                                                                "longitude",
                                                                "sea_water_depth"};

// Register the maker
static RecipeMaker<SeaWaterPotentialTemperature_A> makerSWTempToSWPTemp_(
                    SeaWaterPotentialTemperature_A::Name);

SeaWaterPotentialTemperature_A::SeaWaterPotentialTemperature_A(const Parameters_ & params,
                                    const VaderConfigVars & configVariables):
                                            configVariables_{configVariables}
{
    oops::Log::trace() << "SeaWaterPotentialTemperature_A::SeaWaterPotentialTemperature_A(params)"
        << std::endl;
}

std::string SeaWaterPotentialTemperature_A::name() const
{
    return SeaWaterPotentialTemperature_A::Name;
}

std::string SeaWaterPotentialTemperature_A::product() const
{
    return "sea_water_potential_temperature";
}

std::vector<std::string> SeaWaterPotentialTemperature_A::ingredients() const
{
    return SeaWaterPotentialTemperature_A::Ingredients;
}

size_t SeaWaterPotentialTemperature_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("sea_water_salinity").levels();
}

atlas::FunctionSpace SeaWaterPotentialTemperature_A::productFunctionSpace
                                                (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("sea_water_temperature").functionspace();
}

bool SeaWaterPotentialTemperature_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering SeaWaterPotentialTemperature_A::executeNL function"
        << std::endl;

    // Get fields
    atlas::Field potential_temperature = afieldset.field("sea_water_potential_temperature");
    atlas::Field salinity = afieldset.field("sea_water_salinity");
    atlas::Field latitude = afieldset.field("latitude");
    atlas::Field longitude = afieldset.field("longitude");
    atlas::Field depth = afieldset.field("sea_water_depth");
    atlas::Field insitu_temperature = afieldset.field("sea_water_temperature");

    // Get field views
    auto potential_temperature_view = atlas::array::make_view<double, 2>(potential_temperature);
    auto salinity_view = atlas::array::make_view<double, 2>(salinity);
    auto latitude_view = atlas::array::make_view<double, 2>(latitude);
    auto longitude_view = atlas::array::make_view<double, 2>(longitude);
    auto depth_view = atlas::array::make_view<double, 2>(depth);
    auto insitu_view = atlas::array::make_view<double, 2>(insitu_temperature);

    // Grid dimensions
    size_t grid_size = salinity.shape(0);
    int nlevels = potential_temperature.levels();

    // Local variables
    double pressure, absolute_salinity;

    for (int level = 0; level < nlevels; ++level) {
      for ( size_t jnode = 0; jnode < grid_size ; ++jnode ) {
        // Obtain pressure from depth
        pressure = gsw_p_from_z_f90(-depth_view(jnode, level), latitude_view(jnode, 0));

        // Convert practical salinity to absolute salinity
        absolute_salinity = gsw_sa_from_sp_f90(salinity_view(jnode, level), pressure,
            longitude_view(jnode, 0), latitude_view(jnode, 0));

        // Calculate Sea Water Potential Temperature
        potential_temperature_view(jnode, level) = gsw_pt_from_t_f90(absolute_salinity,
            insitu_view(jnode, level), pressure);
      }
    }

    oops::Log::trace() << "leaving SeaWaterPotentialTemperature_A::executeNL function" << std::endl;

    return true;
}

}  // namespace vader
