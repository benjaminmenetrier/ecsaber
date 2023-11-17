/*
 * (C) Copyright 2022 UCAR.
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "mo/eval_dry_air_density.h"
#include "oops/util/Logger.h"
#include "vader/recipes/DryAirDensityLevelsMinusOne.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char DryAirDensityLevelsMinusOne_B::Name[] = "DryAirDensityLevelsMinusOne_B";
const std::vector<std::string> DryAirDensityLevelsMinusOne_B::
    Ingredients = {"height", "height_levels", "virtual_potential_temperature",
                   "air_pressure_levels"};

// Register the maker
static RecipeMaker<DryAirDensityLevelsMinusOne_B>
    makerDryAirDensityLevelsMinusOne_B_(DryAirDensityLevelsMinusOne_B::Name);

DryAirDensityLevelsMinusOne_B::DryAirDensityLevelsMinusOne_B(const Parameters_ & params,
                                        const VaderConfigVars & configVariables) :
        configVariables_{configVariables}

{
    oops::Log::trace() <<
        "DryAirDensityLevelsMinusOne_B::DryAirDensityLevelsMinusOne_B(params)" << std::endl;
}

std::string DryAirDensityLevelsMinusOne_B::name() const
{
    return DryAirDensityLevelsMinusOne_B::Name;
}

std::string DryAirDensityLevelsMinusOne_B::product() const
{
    return "dry_air_density_levels_minus_one";
}

std::vector<std::string> DryAirDensityLevelsMinusOne_B::ingredients() const
{
    return DryAirDensityLevelsMinusOne_B::Ingredients;
}

size_t DryAirDensityLevelsMinusOne_B::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset["air_pressure_levels"].levels() - 1;
}

atlas::FunctionSpace DryAirDensityLevelsMinusOne_B::productFunctionSpace(const atlas::FieldSet
                                                                 & afieldset) const
{
    return afieldset.field("air_pressure_levels").functionspace();
}

bool DryAirDensityLevelsMinusOne_B::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne_B::executeNL function"
      << std::endl;
    mo::eval_dry_air_density_from_pressure_levels_nl(afieldset);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_B::executeNL function" << std::endl;

    return true;
}

bool DryAirDensityLevelsMinusOne_B::executeTL(atlas::FieldSet & afieldsetTL,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne_B::executeTL function"
        << std::endl;
    mo::eval_dry_air_density_from_pressure_levels_tl(afieldsetTL, afieldsetTraj);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_B::executeTL function" << std::endl;

    return true;
}

bool DryAirDensityLevelsMinusOne_B::executeAD(atlas::FieldSet & afieldsetAD,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne::executeAD function"
        << std::endl;
    mo::eval_dry_air_density_from_pressure_levels_ad(afieldsetAD, afieldsetTraj);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_B::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
