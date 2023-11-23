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
#include "vader/recipes/AirPressureThickness.h"

namespace vader
{

// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPressureThickness_A::Name[] = "AirPressureThickness_A";
const std::vector<std::string> AirPressureThickness_A::Ingredients = {"air_pressure_levels"};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<AirPressureThickness_A>
                    makerAirPressureThickness_A_(AirPressureThickness_A::Name);

// -------------------------------------------------------------------------------------------------

AirPressureThickness_A::AirPressureThickness_A(const AirPressureThickness_AParameters & params,
                                        const VaderConfigVars & configVariables) {
    oops::Log::trace() << "AirPressureThickness_A::AirPressureThickness_A(params)" << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureThickness_A::name() const {
    return AirPressureThickness_A::Name;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureThickness_A::product() const {
    return "air_pressure_thickness";
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> AirPressureThickness_A::ingredients() const {
    return AirPressureThickness_A::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPressureThickness_A::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").levels() - 1;
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressureThickness_A::productFunctionSpace
                                            (const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressureThickness_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering AirPressureThickness_A execute function" << std::endl;
    // Get the fields
    atlas::Field prsi = afieldset.field("air_pressure_levels");
    atlas::Field delp = afieldset.field("air_pressure_thickness");

    // Get the units
    std::string prsi_units, delp_units;
    prsi.metadata().get("units", prsi_units);
    delp.metadata().get("units", delp_units);

    // Assert that the units match
    ASSERT_MSG(prsi_units == delp_units, "In Vader::AirPressureThickness_A::executeNL the "
               "units for pressure " + prsi_units +
               "do not match the pressure thickness units" + delp_units);

    // Set the array views to manipulate the data
    auto prsi_view = atlas::array::make_view<double, 2>(prsi);
    auto delp_view = atlas::array::make_view<double, 2>(delp);

    // Get the grid size
    const int gridSize = delp.shape(0);

    // Compute pressure thickness from pressure at the levels
    for (int level = 0; level < delp.levels(); ++level) {
        for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
            delp_view(jNode, level) = prsi_view(jNode, level+1) - prsi_view(jNode, level);
        }
    }
    oops::Log::trace() << "leaving AirPressureThickness_A execute function" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader
