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
#include "vader/recipes/SurfaceAirPressure.h"

namespace vader
{

// Static attribute initialization
const char SurfaceAirPressure_A::Name[] = "SurfaceAirPressure_A";
const std::vector<std::string> SurfaceAirPressure_A::Ingredients = {"air_pressure_thickness"};

// Register the maker
static RecipeMaker<SurfaceAirPressure_A> makerSurfaceAirPressure_A_(SurfaceAirPressure_A::Name);

// -------------------------------------------------------------------------------------------------

SurfaceAirPressure_A::SurfaceAirPressure_A(const SurfaceAirPressure_AParameters & params,
                                        const VaderConfigVars & configVariables) :
                                        configVariables_{configVariables} {
    oops::Log::trace() << "SurfaceAirPressure_A::SurfaceAirPressure_A Starting" << std::endl;
    oops::Log::trace() << "SurfaceAirPressure_A::SurfaceAirPressure_A Done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string SurfaceAirPressure_A::name() const {
    return SurfaceAirPressure_A::Name;
}

// -------------------------------------------------------------------------------------------------

std::string SurfaceAirPressure_A::product() const {
    return "surface_pressure";
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> SurfaceAirPressure_A::ingredients() const {
    return SurfaceAirPressure_A::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t SurfaceAirPressure_A::productLevels(const atlas::FieldSet & afieldset) const {
    return 1;
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace SurfaceAirPressure_A::productFunctionSpace(const atlas::FieldSet & afieldset)
const {
    return afieldset.field("air_pressure_thickness").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool SurfaceAirPressure_A::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "SurfaceAirPressure_A::executeNL Starting" << std::endl;

    // Get the fields
    atlas::Field delp = afieldset.field("air_pressure_thickness");
    atlas::Field ps = afieldset.field("surface_pressure");

    // Get the units
    std::string delp_units, ps_units, prsi_units;
    delp.metadata().get("units", delp_units);
    ps.metadata().get("units", ps_units);

    // Assert that the units match
    ASSERT_MSG(ps_units == delp_units, "In Vader::SurfaceAirPressure_A::executeNL the units for "
               "surface pressure " + ps_units + "do not match the pressure thickness units"
               + delp_units);

    // Get ptop
    double ptop = configVariables_.getDouble("air_pressure_at_top_of_atmosphere_model");

    // Set the array views to manipulate the data
    auto delp_view = atlas::array::make_view<double, 2>(delp);
    auto ps_view = atlas::array::make_view<double, 2>(ps);

    // Get the grid size
    const int gridSize = delp.shape(0);
    const int nLevel = delp.levels();

    // Set pressure at the surface to surface pressure
    for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
        ps_view(jNode, 0) = ptop;
    }

    // Compute pressure from pressure thickness starting at the surface
    for (int level = 0; level < delp.levels(); ++level) {
        for ( size_t jNode = 0; jNode < gridSize ; ++jNode ) {
            ps_view(jNode, 0) = ps_view(jNode, 0) + delp_view(jNode, level);
        }
    }

    // Return
    oops::Log::trace() << "SurfaceAirPressure_A::executeNL Done" << std::endl;
    return true;
}

}  // namespace vader

// -------------------------------------------------------------------------------------------------
