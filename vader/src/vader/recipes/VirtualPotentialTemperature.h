/*
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_VADER_RECIPES_VIRTUALPOTENTIALTEMPERATURE_H_
#define SRC_VADER_RECIPES_VIRTUALPOTENTIALTEMPERATURE_H_

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

// ------------------------------------------------------------------------------------------------
class VirtualPotentialTemperature_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VirtualPotentialTemperature_AParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

/*! \brief VirtualPotentialTemperature_A class defines a recipe for virtual potential temperature
 *
 *  \details This instantiation of RecipeBase produces virtual potential temperature
 *           using height_levels and hydrostatic exner levels as inputs.
 */
class VirtualPotentialTemperature_A : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef VirtualPotentialTemperature_AParameters Parameters_;

    VirtualPotentialTemperature_A(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override;
    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};

// ------------------------------------------------------------------------------------------------
class VirtualPotentialTemperature_BParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VirtualPotentialTemperature_BParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

/*! \brief VirtualPotentialTemperature_B class defines a recipe for virtual potential temperature
 *
 *  \details This instantiation of RecipeBase produces virtual potential temperature
 *           using specific humidity and potential temperature as inputs.
 */
class VirtualPotentialTemperature_B : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef VirtualPotentialTemperature_BParameters Parameters_;

    VirtualPotentialTemperature_B(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override;
    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};
// ------------------------------------------------------------------------------------------------


}  // namespace vader

#endif  // SRC_VADER_RECIPES_VIRTUALPOTENTIALTEMPERATURE_H_
