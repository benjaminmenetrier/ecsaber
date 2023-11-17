/*
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>
#include <vector>

#include "atlas/field.h"
#include "mo/eval_rain_mixing_ratio.h"
#include "oops/util/Logger.h"
#include "vader/recipes/RainMixingRatio.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char RainMixingRatio_A::Name[] = "RainMixingRatio_A";
const std::vector<std::string> RainMixingRatio_A::
    Ingredients = {"m_r", "m_t"};

// Register the maker
static RecipeMaker<RainMixingRatio_A> makerTotalWater_A_(RainMixingRatio_A::Name);

RainMixingRatio_A::RainMixingRatio_A(const Parameters_ & params,
                                             const VaderConfigVars & configVariables)
{
    oops::Log::trace() <<
        "RainMixingRatio_A::RainMixingRatio_A(params, configVariables)" << std::endl;
}

std::string RainMixingRatio_A::name() const
{
    return RainMixingRatio_A::Name;
}

std::string RainMixingRatio_A::product() const
{
    return "qrain";
}

std::vector<std::string> RainMixingRatio_A::ingredients() const
{
    return RainMixingRatio_A::Ingredients;
}

size_t RainMixingRatio_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return (afieldset["m_t"].levels());
}

atlas::FunctionSpace RainMixingRatio_A::productFunctionSpace(const atlas::FieldSet
                                                             & afieldset) const
{
    return afieldset["m_t"].functionspace();
}

bool RainMixingRatio_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering RainMixingRatio_A::executeNL function"
      << std::endl;
    mo::eval_rain_mixing_ratio_wrt_moist_air_and_condensed_water_nl(afieldset);
    oops::Log::trace() << "leaving RainMixingRatio_A::executeNL function" << std::endl;

    return true;
}

bool RainMixingRatio_A::executeTL(atlas::FieldSet & afieldsetTL,
                                  const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering RainMixingRatio_A::executeTL function"
        << std::endl;
    mo::eval_rain_mixing_ratio_wrt_moist_air_and_condensed_water_tl(afieldsetTL,
                                                                    afieldsetTraj);
    oops::Log::trace() << "leaving RainMixingRatio_A::executeTL function" << std::endl;

    return true;
}

bool RainMixingRatio_A::executeAD(atlas::FieldSet & afieldsetAD,
                                  const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering RainMixingRatio::executeAD function"
        << std::endl;
    mo::eval_rain_mixing_ratio_wrt_moist_air_and_condensed_water_ad(afieldsetAD,
                                                                    afieldsetTraj);
    oops::Log::trace() << "leaving RainMixingRatio_A::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
