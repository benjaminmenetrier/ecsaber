/*
 * (C) Copyright 2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>
#include <map>
#include <vector>

#include "vader/recipe.h"

namespace vader
{

// ------------------------------------------------------------------------------------------------
Recipe::Recipe(const std::string name, const Parameters_ & paramWrapper,
               const VaderConfigVars & configVariables)
    : recipe_{RecipeFactory::create(name, paramWrapper.recipeParams, configVariables)} {}
// ------------------------------------------------------------------------------------------------
std::string Recipe::name() const
{
    return recipe_->name();
}
// ------------------------------------------------------------------------------------------------
std::string Recipe::product() const
{
    return recipe_->product();
}
// ------------------------------------------------------------------------------------------------
std::vector<std::string> Recipe::ingredients() const
{
    return recipe_->ingredients();
}
// ------------------------------------------------------------------------------------------------
size_t Recipe::productLevels(const atlas::FieldSet & afieldset) const
{
    return recipe_->productLevels(afieldset);
}
// ------------------------------------------------------------------------------------------------
atlas::FunctionSpace Recipe::productFunctionSpace(const atlas::FieldSet & afieldset) const
{
    return recipe_->productFunctionSpace(afieldset);
}
// ------------------------------------------------------------------------------------------------
bool Recipe::hasTLAD() const
{
    return recipe_->hasTLAD();
}
// ------------------------------------------------------------------------------------------------
bool Recipe::requiresSetup() const
{
    return recipe_->requiresSetup();
}
// ------------------------------------------------------------------------------------------------
bool Recipe::setup(atlas::FieldSet & afieldset)
{
    return recipe_->setup(afieldset);
}
// ------------------------------------------------------------------------------------------------
bool Recipe::executeNL(atlas::FieldSet & afieldset)
{
    return recipe_->executeNL(afieldset);
}
// ------------------------------------------------------------------------------------------------
bool Recipe::executeTL(atlas::FieldSet & afieldsetTL, const atlas::FieldSet & afieldsetTraj)
{
    return recipe_->executeTL(afieldsetTL, afieldsetTraj);
}
// ------------------------------------------------------------------------------------------------
bool Recipe::executeAD(atlas::FieldSet & afieldsetAD, const atlas::FieldSet & afieldsetTraj)
{
    return recipe_->executeAD(afieldsetAD, afieldsetTraj);
}
// ------------------------------------------------------------------------------------------------
void Recipe::print(std::ostream &os) const
{
    os << *recipe_;
}

}  // namespace vader
