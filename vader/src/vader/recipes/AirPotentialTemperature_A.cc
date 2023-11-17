/*
 * (C) Copyright 2021-2023  UCAR.
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
#include "vader/recipes/AirPotentialTemperature.h"


namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPotentialTemperature_A::Name[] = "AirPotentialTemperature_A";
const std::vector<std::string> AirPotentialTemperature_A::Ingredients = {"air_temperature",
                                                                         "surface_pressure"};

// Register the maker
static RecipeMaker<AirPotentialTemperature_A> makerTempToPTemp_(AirPotentialTemperature_A::Name);

AirPotentialTemperature_A::AirPotentialTemperature_A(const Parameters_ & params,
                                    const VaderConfigVars & configVariables):
                                            configVariables_{configVariables}
{
    oops::Log::trace() << "AirPotentialTemperature_A::AirPotentialTemperature_A(params)"
        << std::endl;
}

std::string AirPotentialTemperature_A::name() const
{
    return AirPotentialTemperature_A::Name;
}

std::string AirPotentialTemperature_A::product() const
{
    // TODO(vahl): update this to air_potential_temperature (models will also need to adapt)
    return "potential_temperature";
}

std::vector<std::string> AirPotentialTemperature_A::ingredients() const
{
    return AirPotentialTemperature_A::Ingredients;
}

size_t AirPotentialTemperature_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("air_temperature").levels();
}

atlas::FunctionSpace AirPotentialTemperature_A::productFunctionSpace
                                                (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("air_temperature").functionspace();
}

bool AirPotentialTemperature_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering AirPotentialTemperature_A::executeNL function"
        << std::endl;

    // Extract values from client config
    const double p0 = configVariables_.getDouble("reference_pressure");
    const double kappa = configVariables_.getDouble("kappa");  // Need better name

    atlas::Field temperature = afieldset.field("air_temperature");
    atlas::Field surface_pressure = afieldset.field("surface_pressure");
    atlas::Field potential_temperature = afieldset.field("potential_temperature");
    std::string t_units, ps_units;

    temperature.metadata().get("units", t_units);
    ASSERT_MSG(t_units == "K", "AirPotentialTemperature_A::executeNL: Incorrect units for "
                            "air_temperature");
    surface_pressure.metadata().get("units", ps_units);
    ASSERT_MSG(ps_units == "Pa", "AirPotentialTemperature_A::executeNL: Incorrect units for "
                            "surface_air_pressure");
    oops::Log::debug() << "AirPotentialTemperature_A::execute: p0 value: " << p0 <<
        std::endl;
    oops::Log::debug() << "AirPotentialTemperature_A::execute: kappa value: " << kappa <<
    std::endl;

    auto temperature_view = atlas::array::make_view<double, 2>(temperature);
    auto surface_pressure_view = atlas::array::make_view<double, 2>(surface_pressure);
    auto potential_temperature_view = atlas::array::make_view<double, 2>(potential_temperature);

    size_t grid_size = surface_pressure.size();

    int nlevels = temperature.levels();
    for (int level = 0; level < nlevels; ++level) {
      for ( size_t jnode = 0; jnode < grid_size ; ++jnode ) {
        potential_temperature_view(jnode, level) =
            temperature_view(jnode, level) * pow(p0 / surface_pressure_view(jnode, 0), kappa);
      }
    }

    oops::Log::trace() << "leaving AirPotentialTemperature_A::executeNL function" << std::endl;

    return true;
}

}  // namespace vader
