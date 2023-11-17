/*
 * (C) Copyright 2021-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_VADER_RECIPES_VWIND_AT_10M_H_
#define SRC_VADER_RECIPES_VWIND_AT_10M_H_

#include <map>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

class vwind_at_10m_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(vwind_at_10m_AParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

// ------------------------------------------------------------------------------------------------
/*! \brief vwind_at_10m_A class defines a recipe for surface northward wind at 10 m
 *
 *  \details This instantiation of RecipeBase produces dry surface northward wind at 10 m
 *           using northward wind at the lowest atmospheric level as input.
 */
class vwind_at_10m_A : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef vwind_at_10m_AParameters Parameters_;

    vwind_at_10m_A(const Parameters_ &, const VaderConfigVars &);

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
};

}  // namespace vader

#endif  // SRC_VADER_RECIPES_VWIND_AT_10M_H_
