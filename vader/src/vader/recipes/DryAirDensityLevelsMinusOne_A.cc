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
const char DryAirDensityLevelsMinusOne_A::Name[] = "DryAirDensityLevelsMinusOne_A";
const std::vector<std::string> DryAirDensityLevelsMinusOne_A::
    Ingredients = {"height", "height_levels", "virtual_potential_temperature",
                   "air_pressure_levels_minus_one"};

// Register the maker
static RecipeMaker<DryAirDensityLevelsMinusOne_A>
    makerDryAirDensityLevelsMinusOne_A_(DryAirDensityLevelsMinusOne_A::Name);

DryAirDensityLevelsMinusOne_A::DryAirDensityLevelsMinusOne_A(const Parameters_ & params,
                                        const VaderConfigVars & configVariables) :
        configVariables_{configVariables}

{
    oops::Log::trace() <<
        "DryAirDensityLevelsMinusOne_A::DryAirDensityLevelsMinusOne_A(params)" << std::endl;
}

std::string DryAirDensityLevelsMinusOne_A::name() const
{
    return DryAirDensityLevelsMinusOne_A::Name;
}

std::string DryAirDensityLevelsMinusOne_A::product() const
{
    return "dry_air_density_levels_minus_one";
}

std::vector<std::string> DryAirDensityLevelsMinusOne_A::ingredients() const
{
    return DryAirDensityLevelsMinusOne_A::Ingredients;
}

size_t DryAirDensityLevelsMinusOne_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset["air_pressure_levels_minus_one"].levels();
}

atlas::FunctionSpace DryAirDensityLevelsMinusOne_A::productFunctionSpace(const atlas::FieldSet
                                                                 & afieldset) const
{
    return afieldset.field("air_pressure_levels_minus_one").functionspace();
}

bool DryAirDensityLevelsMinusOne_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne_A::executeNL function"
      << std::endl;
    mo::eval_dry_air_density_from_pressure_levels_minus_one_nl(afieldset);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_A::executeNL function" << std::endl;

    return true;
}

bool DryAirDensityLevelsMinusOne_A::executeTL(atlas::FieldSet & afieldsetTL,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne_A::executeTL function"
        << std::endl;
    mo::eval_dry_air_density_from_pressure_levels_minus_one_tl(afieldsetTL, afieldsetTraj);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_A::executeTL function" << std::endl;

    return true;
}

bool DryAirDensityLevelsMinusOne_A::executeAD(atlas::FieldSet & afieldsetAD,
                                      const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering DryAirDensityLevelsMinusOne::executeAD function"
        << std::endl;
    mo::eval_dry_air_density_from_pressure_levels_minus_one_ad(afieldsetAD, afieldsetTraj);
    oops::Log::trace() << "leaving DryAirDensityLevelsMinusOne_A::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
