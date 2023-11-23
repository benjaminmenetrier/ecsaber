/*
 * (C) Crown Copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_VADER_RECIPES_HYDROSTATICEXNERLEVELS_H_
#define SRC_VADER_RECIPES_HYDROSTATICEXNERLEVELS_H_

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

class HydrostaticExnerLevels_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerLevels_AParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

// ------------------------------------------------------------------------------------------------
/*! \brief HydrostaticExnerLevels_A class defines a recipe for hydrostatic exner levels
 *
 *  \details This instantiation of RecipeBase produces hydrostatic exner levels
 *           using height_levels, virtual_potential_temperature and
 *           air_pressure_levels_minus_one as inputs.
 */
class HydrostaticExnerLevels_A : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef HydrostaticExnerLevels_AParameters Parameters_;

    HydrostaticExnerLevels_A(const Parameters_ &, const VaderConfigVars &);

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

}  // namespace vader

#endif  // SRC_VADER_RECIPES_HYDROSTATICEXNERLEVELS_H_
