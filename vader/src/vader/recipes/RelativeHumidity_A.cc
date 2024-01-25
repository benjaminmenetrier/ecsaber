/*
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "mo/eval_relative_humidity.h"
#include "oops/util/Logger.h"
#include "vader/recipes/RelativeHumidity.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char RelativeHumidity_A::Name[] = "RelativeHumidity_A";
const std::vector<std::string> RelativeHumidity_A::Ingredients = {"air_temperature",
                                                                  "dlsvpdT",
                                                                  "specific_humidity",
                                                                  "qsat"};

// Register the maker
static RecipeMaker<RelativeHumidity_A> makerRelativeHumidity_(RelativeHumidity_A::Name);

RelativeHumidity_A::RelativeHumidity_A(const Parameters_ & params,
                                       const VaderConfigVars & configVariables) :
    configVariables_{configVariables}
{
    oops::Log::trace() << "RelativeHumidity_A::RelativeHumidity_A(params)" << std::endl;
}

std::string RelativeHumidity_A::name() const
{
    return RelativeHumidity_A::Name;
}

std::string RelativeHumidity_A::product() const
{
    return "relative_humidity";
}

std::vector<std::string> RelativeHumidity_A::ingredients() const
{
    return RelativeHumidity_A::Ingredients;
}

size_t RelativeHumidity_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("air_temperature").levels();
}

atlas::FunctionSpace RelativeHumidity_A::productFunctionSpace
                                              (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("air_temperature").functionspace();
}

bool RelativeHumidity_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering RelativeHumidity_A::executeNL function" << std::endl;

    // Actually requiring specific_humidity and qsat only
    const bool rvalue = mo::eval_relative_humidity_nl(afieldset);

    oops::Log::trace() << "leaving RelativeHumidity_A::executeNL function" << std::endl;

    return rvalue;
}

bool RelativeHumidity_A::executeTL(atlas::FieldSet & afieldsetTL,
                                   const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering RelativeHumidity_A::executeTL function" << std::endl;

    // Requiring air_temperature and specific_humidity in perturbation fieldSet,
    // But specific_humidity, qsat and dlsvpdt in trajectory fieldSet.
    mo::eval_relative_humidity_tl(afieldsetTL, afieldsetTraj);

    oops::Log::trace() << "leaving RelativeHumidity_A::executeTL function" << std::endl;

    return true;
}

bool RelativeHumidity_A::executeAD(atlas::FieldSet & afieldsetAD,
                                   const atlas::FieldSet & afieldsetTraj)
{
    oops::Log::trace() << "entering RelativeHumidity_A::executeAD function" << std::endl;

    mo::eval_relative_humidity_ad(afieldsetAD, afieldsetTraj);

    oops::Log::trace() << "leaving RelativeHumidity_A::executeAD function" << std::endl;

    return true;
}

}  // namespace vader
