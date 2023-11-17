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
#include "vader/recipes/LnAirPressureAtInterface.h"

namespace vader {


// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char LnAirPressureAtInterface_A::Name[] = "LnAirPressureAtInterface_A";
const std::vector<std::string> LnAirPressureAtInterface_A::Ingredients = {"air_pressure_levels"};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<LnAirPressureAtInterface_A>
       makerLnAirPressureAtInterface_A_(LnAirPressureAtInterface_A::Name);

// -------------------------------------------------------------------------------------------------

LnAirPressureAtInterface_A::LnAirPressureAtInterface_A(const
                                                     LnAirPressureAtInterface_AParameters & params,
                                                     const VaderConfigVars & configVariables)
{
    oops::Log::trace() << "LnAirPressureAtInterface_A::LnAirPressureAtInterface_A Starting"
                       << std::endl;
    oops::Log::trace() << "LnAirPressureAtInterface_A::LnAirPressureAtInterface_A Done"
                       << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string LnAirPressureAtInterface_A::name() const {
    return LnAirPressureAtInterface_A::Name;
}

// -------------------------------------------------------------------------------------------------

std::string LnAirPressureAtInterface_A::product() const {
    return "ln_air_pressure_at_interface";
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> LnAirPressureAtInterface_A::ingredients() const {
    return LnAirPressureAtInterface_A::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t LnAirPressureAtInterface_A::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").levels();
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace LnAirPressureAtInterface_A::productFunctionSpace(const
                                                                atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool LnAirPressureAtInterface_A::executeNL(atlas::FieldSet & afieldset) {
    //
    oops::Log::trace() << "LnAirPressureAtInterface_A::executeNL Starting" << std::endl;

    // Get fields
    atlas::Field airPressureLevelsF = afieldset.field("air_pressure_levels");
    atlas::Field lnAirPressureAtInterfaceF = afieldset.field("ln_air_pressure_at_interface");

    // Get field views
    auto airPressureLevels = atlas::array::make_view<double, 2>(airPressureLevelsF);
    auto lnAirPressureAtInterface = atlas::array::make_view<double, 2>(lnAirPressureAtInterfaceF);

    // Grid dimensions
    size_t h_size = airPressureLevelsF.shape(0);
    int v_size = airPressureLevelsF.levels();

    // Calculate the output variable
    for (int vv = 0; vv < v_size; ++vv) {
      for ( size_t hh = 0; hh < h_size ; ++hh ) {
        lnAirPressureAtInterface(hh, vv) = log(airPressureLevels(hh, vv));
      }
    }

    // Return
    oops::Log::trace() << "LnAirPressureAtInterface_A::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
