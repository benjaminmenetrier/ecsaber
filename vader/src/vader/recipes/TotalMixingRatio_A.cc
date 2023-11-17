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
#include "mo/eval_total_mixing_ratio.h"
#include "oops/util/Logger.h"
#include "vader/recipes/TotalMixingRatio.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char TotalMixingRatio_A::Name[] = "TotalMixingRatio_A";
const std::vector<std::string> TotalMixingRatio_A::
    Ingredients = {"m_v",   // CCPP name: water_vapor_mixing_ratio_wrt_dry_air
                   "m_ci",  // CCPP name: cloud_ice_mixing_ratio_wrt_dry_air
                   "m_cl",  // CCPP name: cloud_liquid_water_mixing_ratio_wrt_dry_air
                   "m_r"};  // CCPP name: rain_mixing_ratio_wrt_dry_air

// Register the maker
static RecipeMaker<TotalMixingRatio_A> makerTotalMixingRatio_A_(TotalMixingRatio_A::Name);

TotalMixingRatio_A::TotalMixingRatio_A(const Parameters_ & params,
                           const VaderConfigVars & configVariables)
{
    oops::Log::trace() <<
        "TotalMixingRatio_A::TotalMixingRatio_A(params, configVariables)" << std::endl;
}

std::string TotalMixingRatio_A::name() const
{
    return TotalMixingRatio_A::Name;
}

std::string TotalMixingRatio_A::product() const
{
    return "m_t";  // CCPP name: total_water_mixing_ratio_wrt_dry_air
}

std::vector<std::string> TotalMixingRatio_A::ingredients() const
{
    return TotalMixingRatio_A::Ingredients;
}

size_t TotalMixingRatio_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return (afieldset["m_v"].levels());
}

atlas::FunctionSpace TotalMixingRatio_A::productFunctionSpace(const atlas::FieldSet
                                                        & afieldset) const
{
    return afieldset["m_v"].functionspace();
}

bool TotalMixingRatio_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering TotalMixingRatio_A::executeNL function"
      << std::endl;
    mo::eval_total_mixing_ratio_nl(afieldset);
    oops::Log::trace() << "leaving TotalMixingRatio_A::executeNL function" << std::endl;

    return true;
}

bool TotalMixingRatio_A::executeTL(atlas::FieldSet & afieldsetTL,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering TotalMixingRatio_A::executeTL function"
        << std::endl;
    mo::eval_total_mixing_ratio_tl(afieldsetTL, afieldsetTraj);
    oops::Log::trace() << "leaving TotalMixingRatio_A::executeTL function" << std::endl;

    return true;
}

bool TotalMixingRatio_A::executeAD(atlas::FieldSet & afieldsetAD,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering TotalMixingRatio::executeAD function"
        << std::endl;
    mo::eval_total_mixing_ratio_ad(afieldsetAD, afieldsetTraj);
    oops::Log::trace() << "leaving TotalMixingRatio_A::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
