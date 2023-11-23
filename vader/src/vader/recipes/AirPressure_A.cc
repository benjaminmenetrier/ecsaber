/*
 * (C) Copyright 2023  UCAR.
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
#include "vader/recipes/AirPressure.h"

namespace vader {

// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPressure_A::Name[] = "AirPressure_A";
const std::vector<std::string> AirPressure_A::Ingredients = {"air_pressure_levels"};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<AirPressure_A> makerAirPressure_A_(AirPressure_A::Name);

// -------------------------------------------------------------------------------------------------

AirPressure_A::AirPressure_A(const AirPressure_AParameters & params,
                             const VaderConfigVars & configVariables) :
                                    configVariables_{configVariables}
{
    oops::Log::trace() << "AirPressure_A::AirPressure_A Starting" << std::endl;
    oops::Log::trace() << "AirPressure_A::AirPressure_A Done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressure_A::name() const {
    return AirPressure_A::Name;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressure_A::product() const {
    return "air_pressure";
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> AirPressure_A::ingredients() const {
    return AirPressure_A::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPressure_A::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").levels() - 1;
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressure_A::productFunctionSpace(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressure_A::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "AirPressure_A::executeNL Starting" << std::endl;

    // Extract value from client config
    const double kappa = configVariables_.getDouble("kappa");  // Need better name

    // Get fields
    atlas::Field airPressureLevelsF = afieldset.field("air_pressure_levels");
    atlas::Field airPressureF = afieldset.field("air_pressure");

    auto airPressureLevels = atlas::array::make_view<double, 2>(airPressureLevelsF);
    auto airPressure = atlas::array::make_view<double, 2>(airPressureF);

    // Grid dimensions
    size_t h_size = airPressureLevelsF.shape(0);
    int v_size = airPressureLevelsF.levels() - 1;

    // kappa variations
    const double kap1 = kappa + 1.0;
    const double kapr = 1.0 / kappa;

    // Calculate the output variable
    for (int vv = 0; vv < v_size; ++vv) {
      for ( size_t hh = 0; hh < h_size ; ++hh ) {
        airPressure(hh, vv) = pow(((pow(airPressureLevels(hh, vv+1), kap1) -
                                    pow(airPressureLevels(hh, vv), kap1)) /
                                    (kap1*(airPressureLevels(hh, vv+1) -
                                    airPressureLevels(hh, vv)))), kapr);
      }
    }

    // Return
    oops::Log::trace() << "AirPressure_A::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
