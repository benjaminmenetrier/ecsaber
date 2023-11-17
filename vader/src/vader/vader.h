/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"
#include "RecipeBase.h"
#include "VaderParameters.h"

namespace vader {

typedef std::map<std::string, std::vector<std::string>> cookbookConfigType;
// configCookbookKey is a key that can optionally be provided in the eckit::LocalConfiguration
// that can be passed as the second parameter of the Vader constructor. If passed, the value
// is used to define the Vader cookbook for the vader instance. If not passed, the default Vader
// cookbook defined in DefaultCookbook.h will be used.
static const char configCookbookKey[]  = "cookbook";

// configModelVarsKey is a key that can optionally be provided in the eckit::LocalConfiguration
// that can be passed as the second parameter of the Vader constructor. If passed, the values
// will be stored in Vader's configVariables_ property. Whenever a recipe is executed that
// requires a model-specific value in its algorithm, the recipe will look in that Configuration for
// the value, and should generate an error if it isn't found. (Instead of using a default value.)
static const char configModelVarsKey[] = "model data";

// ------------------------------------------------------------------------------------------------
/*! \brief Vader class to handle variable transformations
 *
 *  \details This class provides generic variable transformations via the
 *           changeVar method which is passed an atlas fieldset.
 *
 *           Throughout the Vader code, the primary variables are named using a
 *           metaphor involving ingredients, recipies, and a cookbook.
 *
 *           A 'recipe' is is an object that can produce a single output
 *           variable when provided with a list of required input variables.
 *           The input variables are referred to as the 'ingredients' to the
 *           recipe. The 'cookbook' is the container (a map) that
 *           contains the recipes to be attempted when specified output variable
 *           is desired. The cookbook can contain multiple recipes that produce
 *           the same output variable.
 */

class Vader  : public util::Printable {
 public:
    static const std::string classname() {return "Vader";}
    static const cookbookConfigType defaultCookbookDefinition;
    typedef  std::vector<std::pair<std::string,
                                   const std::unique_ptr<RecipeBase> & >> vaderPlanType;
    Vader(const VaderParameters & parameters,
          const eckit::Configuration & config = eckit::LocalConfiguration());
    Vader(const Vader &) = delete;
    Vader& operator=(const Vader &) = delete;
    ~Vader();

    std::vector<std::string> getPlanNames() const;
    std::vector<std::string> getPlanNames(vaderPlanType plan) const;

    /// Calculates as many variables in the list as possible
    oops::patch::Variables changeVar(atlas::FieldSet &, oops::patch::Variables &,
                              vaderPlanType plan = vaderPlanType()) const;
    oops::patch::Variables changeVarTraj(atlas::FieldSet &, oops::patch::Variables &);
    oops::patch::Variables changeVarTL(atlas::FieldSet &, oops::patch::Variables &) const;
    oops::patch::Variables changeVarAD(atlas::FieldSet &, oops::patch::Variables &) const;

 private:
    std::map<std::string, std::vector<std::unique_ptr<RecipeBase>>> cookbook_;
    vaderPlanType recipeExecutionPlan_;
    atlas::FieldSet trajectory_;
    const eckit::LocalConfiguration configVariables_;
    std::map<std::string, std::vector<std::string>>
        getDefaultCookbookDef();

    void createCookbook(const cookbookConfigType &,
                        const std::vector<RecipeParametersWrapper> & allRecpParamWraps =
                              std::vector<RecipeParametersWrapper>());

    bool planVariable(std::vector<std::string> &,
                      oops::patch::Variables &,
                      const std::string,
                      bool,
                      oops::patch::Variables &,
                      vaderPlanType &) const;

    void executePlanNL(atlas::FieldSet &, const vaderPlanType &) const;
    void executePlanTL(atlas::FieldSet &, const vaderPlanType &) const;
    void executePlanAD(atlas::FieldSet &, const vaderPlanType &) const;
    void print(std::ostream &) const;
};

}  // namespace vader
