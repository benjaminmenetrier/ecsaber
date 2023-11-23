 /*
 * (C) Copyright 2021-2022  UCAR.
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
#include "vader/recipes/AirPressureAtInterface.h"

namespace vader
{

// Static attribute initialization
const char AirPressureAtInterface_B::Name[] = "AirPressureAtInterface_B";
const std::vector<std::string> AirPressureAtInterface_B::Ingredients = {"air_pressure_thickness"};

// Register the maker
static RecipeMaker<AirPressureAtInterface_B>
       makerAirPressureAtInterface_B_(AirPressureAtInterface_B::Name);

// -------------------------------------------------------------------------------------------------

AirPressureAtInterface_B::AirPressureAtInterface_B(
                                                 const AirPressureAtInterface_BParameters & params,
                                        const VaderConfigVars & configVariables) :
                                        configVariables_{configVariables}
{
    oops::Log::trace() << "AirPressureAtInterface_B::AirPressureAtInterface_B(params) Starting"
                       << std::endl;
    oops::Log::trace() << "AirPressureAtInterface_B::AirPressureAtInterface_B(params) Done"
                       << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureAtInterface_B::name() const {
    return AirPressureAtInterface_B::Name;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureAtInterface_B::product() const {
    return "air_pressure_levels";
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> AirPressureAtInterface_B::ingredients() const {
    return AirPressureAtInterface_B::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPressureAtInterface_B::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_thickness").levels()+1;
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressureAtInterface_B::productFunctionSpace
                                                        (const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_thickness").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressureAtInterface_B::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "AirPressureAtInterface_B::executeNL Starting" << std::endl;

    double ptop = configVariables_.getDouble("air_pressure_at_top_of_atmosphere_model");

    // Get the fields
    atlas::Field delp = afieldset.field("air_pressure_thickness");
    atlas::Field prsi = afieldset.field("air_pressure_levels");

    // Get the units
    std::string delp_units, prsi_units;
    delp.metadata().get("units", delp_units);
    prsi.metadata().get("units", prsi_units);

    // Assert that the units match
    ASSERT_MSG(prsi_units == delp_units, "In Vader::AirPressureAtInterface_B::executeNL the units "
               "for pressure at the levels" + prsi_units + "do not match the pressure thickness "
               "units " + delp_units);

    // Set the array views to manipulate the data
    auto delp_view = atlas::array::make_view<double, 2>(delp);
    auto prsi_view = atlas::array::make_view<double, 2>(prsi);

    // Get the grid size
    const int gridSize = delp.shape(0);
    const int nLevel = prsi.levels() - 1;  // Reduce by 1 since index begins at 0

    // Set pressure at the surface to surface pressure
    for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
        prsi_view(jNode, 0) = ptop;
    }

    // Compute pressure from pressure thickness starting at the surface
    for (int level = 1; level < prsi.levels(); ++level) {
        for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
            prsi_view(jNode, level) = prsi_view(jNode, level - 1) + delp_view(jNode, level - 1);
        }
    }

    // Return
    oops::Log::trace() << "AirPressureAtInterface_B::executeNL Done" << std::endl;
    return true;
}

}  // namespace vader

// -------------------------------------------------------------------------------------------------
