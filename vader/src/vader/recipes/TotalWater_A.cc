/*
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>
#include <vector>

#include "atlas/field.h"
// #include "atlas/util/Metadata.h"
#include "mo/eval_moisture_incrementing_operator.h"
#include "oops/util/Logger.h"
#include "vader/recipes/TotalWater.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char TotalWater_A::Name[] = "TotalWater_A";
const std::vector<std::string> TotalWater_A::
    Ingredients = {"specific_humidity",
                   "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
                   "mass_content_of_cloud_ice_in_atmosphere_layer"};

// Register the maker
static RecipeMaker<TotalWater_A> makerTotalWater_A_(TotalWater_A::Name);

TotalWater_A::TotalWater_A(const Parameters_ & params,
                           const VaderConfigVars & configVariables)
{
    oops::Log::trace() <<
        "TotalWater_A::TotalWater_A(params, configVariables)" << std::endl;
}

std::string TotalWater_A::name() const
{
    return TotalWater_A::Name;
}

std::string TotalWater_A::product() const
{
    return "qt";
}

std::vector<std::string> TotalWater_A::ingredients() const
{
    return TotalWater_A::Ingredients;
}

size_t TotalWater_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return (afieldset["specific_humidity"].levels());
}

atlas::FunctionSpace TotalWater_A::productFunctionSpace(const atlas::FieldSet
                                                        & afieldset) const
{
    return afieldset["specific_humidity"].functionspace();
}

bool TotalWater_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering TotalWater_A::executeNL function"
      << std::endl;
    // There is no non-linear version of this, the full version
    // is given by the tangent linear.
    mo::eval_total_water_tl(afieldset, afieldset);
    oops::Log::trace() << "leaving TotalWater_A::executeNL function" << std::endl;

    return true;
}

bool TotalWater_A::executeTL(atlas::FieldSet & afieldsetTL,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering TotalWater_A::executeTL function"
        << std::endl;
    mo::eval_total_water_tl(afieldsetTL, afieldsetTraj);
    oops::Log::trace() << "leaving TotalWater_A::executeTL function" << std::endl;

    return true;
}

bool TotalWater_A::executeAD(atlas::FieldSet & afieldsetAD,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering TotalWater::executeAD function"
        << std::endl;
    mo::eval_total_water_ad(afieldsetAD, afieldsetTraj);
    oops::Log::trace() << "leaving TotalWater_A::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
