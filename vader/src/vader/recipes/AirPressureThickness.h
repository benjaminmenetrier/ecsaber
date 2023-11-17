/*
 * (C) Copyright 2021-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "vader/RecipeBase.h"

namespace vader
{

// -------------------------------------------------------------------------------------------------

class AirPressureThickness_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirPressureThickness_AParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
     // No parameters for this recipe currently
};

class AirPressureThickness_A : public RecipeBase
{
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef AirPressureThickness_AParameters Parameters_;

    AirPressureThickness_A(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override;
    std::string product() const override;
    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader

