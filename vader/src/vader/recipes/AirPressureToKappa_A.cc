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
#include "vader/recipes/AirPressureToKappa.h"

namespace vader {

// -------------------------------------------------------------------------------------------------

// Static attribute initialization
const char AirPressureToKappa_A::Name[] = "AirPressureToKappa_A";
const std::vector<std::string> AirPressureToKappa_A::Ingredients = {"air_pressure_levels",
                                                                    "ln_air_pressure_at_interface"};

// -------------------------------------------------------------------------------------------------

// Register the maker
static RecipeMaker<AirPressureToKappa_A> makerAirPressureToKappa_A_(AirPressureToKappa_A::Name);

// -------------------------------------------------------------------------------------------------

AirPressureToKappa_A::AirPressureToKappa_A(const AirPressureToKappa_AParameters & params,
                                        const VaderConfigVars & configVariables)
: configVariables_{configVariables} {
    oops::Log::trace() << "AirPressureToKappa_A::AirPressureToKappa_A Starting" << std::endl;
    oops::Log::trace() << "AirPressureToKappa_A::AirPressureToKappa_A Done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureToKappa_A::name() const {
    return AirPressureToKappa_A::Name;
}

// -------------------------------------------------------------------------------------------------

std::string AirPressureToKappa_A::product() const {
    return "air_pressure_to_kappa";
}

// -------------------------------------------------------------------------------------------------

std::vector<std::string> AirPressureToKappa_A::ingredients() const {
    return AirPressureToKappa_A::Ingredients;
}

// -------------------------------------------------------------------------------------------------

size_t AirPressureToKappa_A::productLevels(const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").levels() - 1;
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace AirPressureToKappa_A::productFunctionSpace(
  const atlas::FieldSet & afieldset) const {
    return afieldset.field("air_pressure_levels").functionspace();
}

// -------------------------------------------------------------------------------------------------

bool AirPressureToKappa_A::executeNL(atlas::FieldSet & afieldset) {
    oops::Log::trace() << "AirPressureToKappa_A::executeNL Starting" << std::endl;

    const double kappa = configVariables_.getDouble("kappa");  // Need better name

    // Get fields
    atlas::Field airPressureLevelsF = afieldset.field("air_pressure_levels");
    atlas::Field lnAirPressureAtInterfaceF = afieldset.field("ln_air_pressure_at_interface");
    atlas::Field airPressureToKappaF = afieldset.field("air_pressure_to_kappa");

    // Get field views
    auto airPressureLevels = atlas::array::make_view<double, 2>(airPressureLevelsF);
    auto lnAirPressureAtInterface = atlas::array::make_view<double, 2>(lnAirPressureAtInterfaceF);
    auto airPressureToKappa = atlas::array::make_view<double, 2>(airPressureToKappaF);

    // Grid dimensions
    size_t h_size = airPressureLevelsF.shape(0);
    int v_size = airPressureLevelsF.levels() - 1;

    // Local variables
    double pk1, pk2;

    // Calculate the output variable
    for (int vv = 0; vv < v_size; ++vv) {
      for ( size_t hh = 0; hh < h_size ; ++hh ) {
        // Temporary exp(kappa(ln(p)))
        pk1 = exp(kappa*lnAirPressureAtInterface(hh, vv));
        pk2 = exp(kappa*lnAirPressureAtInterface(hh, vv+1));
        // Compute p to the kappa
        airPressureToKappa(hh, vv) = (pk2 - pk1) /
                 (kappa*(lnAirPressureAtInterface(hh, vv+1) - lnAirPressureAtInterface(hh, vv)));
      }
    }

    // Return
    oops::Log::trace() << "AirPressureToKappa_A::executeNL Done" << std::endl;
    return true;
}

// -------------------------------------------------------------------------------------------------

}  // namespace vader

// -------------------------------------------------------------------------------------------------
