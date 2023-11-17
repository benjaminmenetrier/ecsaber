/*
 * (C) Copyright 2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/option/Options.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"
#include "vader/DefaultCookbook.h"
#include "vader/vader.h"

namespace vader {

// ------------------------------------------------------------------------------------------------
Vader::~Vader() {
    oops::Log::trace() << "Vader::~Vader done" << std::endl;
}
// ------------------------------------------------------------------------------------------------
void Vader::createCookbook(const cookbookConfigType & cookbook,
                           const std::vector<RecipeParametersWrapper> & allRecpParamWraps) {
    oops::Log::trace() << "entering Vader::createCookbook" << std::endl;
    std::vector<std::unique_ptr<RecipeBase>> recipes;
    for (auto defEntry : cookbook) {
        recipes.clear();
        for (auto recipeName : defEntry.second) {
            // There might not be any recipe parameters at all.
            // There might not be parameters for THIS recipe.
            // We must prepare for all eventualities.
            bool parametersFound = false;
            for (auto & singleRecpParamWrap : allRecpParamWraps) {
                if (singleRecpParamWrap.recipeParams.value().name.value() == recipeName) {
                    recipes.push_back(std::unique_ptr<RecipeBase>
                        (RecipeFactory::create(recipeName,
                                               singleRecpParamWrap.recipeParams,
                                               configVariables_)));
                    parametersFound = true;
                    break;
                }
            }
            if (!parametersFound) {
                auto emptyRecipeParams = RecipeFactory::createParameters(recipeName);
                recipes.push_back(std::unique_ptr<RecipeBase>
                                  (RecipeFactory::create(recipeName, *emptyRecipeParams,
                                                         configVariables_)));
            }
        }
        cookbook_[defEntry.first] = std::move(recipes);
    }
    oops::Log::trace() << "leaving Vader::createCookbook" << std::endl;
}
// ------------------------------------------------------------------------------------------------
Vader::Vader(const VaderParameters & parameters, const eckit::Configuration & config) :
             configVariables_(config.getSubConfiguration(configModelVarsKey)) {
    util::Timer timer(classname(), "Vader");
    oops::Log::trace() << "entering Vader::Vader(parameters, constructConfig) " << std::endl;

    // Vader is designed to function without parameters. So VaderParameters
    // should not have any RequiredParameters.
    //
    // To simplify things for vader clients, they should declare vader Parameters with a
    // default construction of empty/default VaderParameters. i.e. their Parameters should contain:
    // oops::Parameter<vader::VaderParameters> vader{"vader", {}, this};
    //
    cookbookConfigType cookbookDefinition;
    if (config.has(configCookbookKey)) {
        // Convert eckit::Configuration to cookbookConfigType (map)
        auto configCookbookDefinition = config.getSubConfiguration(configCookbookKey);
        for (auto configCookbookKey : configCookbookDefinition.keys()) {
            cookbookDefinition[configCookbookKey] =
                                        configCookbookDefinition.getStringVector(configCookbookKey);
        }
    } else {
        cookbookDefinition = Vader::defaultCookbookDefinition;
    }
    if (parameters.recipeParams.value() == boost::none) {
        createCookbook(cookbookDefinition);
    } else {
        createCookbook(cookbookDefinition, *parameters.recipeParams.value());
    }
    oops::Log::trace() << "leaving Vader::Vader(parameters, constructConfig) " << std::endl;
}
// ------------------------------------------------------------------------------------------------
std::vector<std::string> Vader::getPlanNames() const {
    return getPlanNames(recipeExecutionPlan_);
}
// ------------------------------------------------------------------------------------------------
/*! \brief Get Plan Names
*
* \details **getPlanNames** provides a list of the names, in order, of the recipes in
* plan. (This is primarily useful for testing the Vader algorithm.)
*
* \param[in] plan A populated variable of type vaderPlanType
* \returns Ordered list of the names of the recipes in plan.
*/
std::vector<std::string> Vader::getPlanNames(const vaderPlanType plan) const {
    std::vector<std::string> recipeNames;
    for (const auto & varPlan : plan)
    {
        recipeNames.push_back(varPlan.second->name());
    }
    return recipeNames;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Change Variable
*
* \details **changeVar** is called externally to invoke Vader's non-linear variable change
* functionality. The caller passes an Atlas FieldSet that contains only fields that have already
* been populated with values. These fields serve as the initial ingredients for recipes. The names
* of the desired output variables are passed via the neededVars parameter. After this method is
* complete, Vader will have popluated all the desired variables it can based on
* the ingredients it was given and the recipes in its cookbook. The names of the
* variables Vader created are returned to the caller in two different ways:
* * they are removed from neededVars
* * they are put into the method return object.
* Note that Vader may create some variables not listed in neededVars if they are required
* intermediate ingredients for recipes in order to produce one or more of the neededVars.
* Any variable names remaining in neededVars could not be created by Vader.
*
* \param[in,out] afieldset This is the FieldSet described above
* \param[in,out] neededVars Names of variables to try to produce
* \param[out] plan Optional arg allowing the caller to recieve the recipe execution plan used
                   (can be used for testing the planVariable routine)
* \returns List of variables VADER was able to produce
*
*/
oops::patch::Variables Vader::changeVar(atlas::FieldSet & afieldset,
                                 oops::patch::Variables & neededVars,
                                 vaderPlanType plan) const {
    util::Timer timer(classname(), "changeVar");
    oops::Log::trace() << "entering Vader::changeVar " << std::endl;
    oops::Log::info() << "Variables requested from Vader (non-linear): " << neededVars.variables()
        << std::endl;

    // onlyIngredientVars are the variables that Vader should NOT later return in varsProduced
    // (i.e. variables both in afieldset AND neededVars will be returned in varsProduced)
    oops::patch::Variables onlyIngredientVars(afieldset.field_names());
    onlyIngredientVars -= neededVars;
    plan.clear();

    auto ingredientNames = afieldset.field_names();
    oops::Log::info() << "Variables provided to Vader: " << ingredientNames << std::endl;
    // Loop through all the requested fields in neededVars
    // Since neededVars can be modified by planVariable and planVariable calls
    // itself recursively, we make a copy of the list here before we start.
    std::vector<std::string> targetVariables{neededVars.variables()};
    bool recipesNeedTLAD = false;  // It's OK here to plan recipes with no TL/AD methods

    for (auto targetVariable : targetVariables) {
        oops::Log::debug() <<
            "Vader::changeVar processing top-level requested variable: "
            << targetVariable << std::endl;
        // Since this function is called recursively, make sure targetVariable is
        // still needed
        if (!neededVars.has(targetVariable)) {
            oops::Log::debug() << targetVariable <<
                " is no longer in the variable list neededVars. Moving on." << std::endl;
            continue;
        }
        oops::patch::Variables excludedVars;
        planVariable(ingredientNames, neededVars, targetVariable, recipesNeedTLAD,
                     excludedVars, plan);
    }
    executePlanNL(afieldset, plan);

    if (neededVars.variables().size()) {
        oops::Log::info() << "Requested variables Vader could not produce: "
            << neededVars.variables() << std::endl;
    } else {
        oops::Log::info() << "Vader produced all requested variables."  << std::endl;
    }

    oops::patch::Variables varsProduced(afieldset.field_names());
    varsProduced -= onlyIngredientVars;
    oops::Log::trace() << "leaving Vader::changeVar" << std::endl;
    return varsProduced;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Change Variable-Trajectory
*
* \details **changeVarTraj** is called externally to set up the trajectory for subsequent calls to
* vader's changeVarTL and changeVarAD methods. (The corresponding method used to be called
* 'setTrajectory' in OOPS.) It performs the same non-linear variable change logic
* as the changeVar method, but also saves the result in Vader's trajectory_ property.
*
* \param[in,out] afieldset This is the input/output FieldSet, same as in changeVar.
* \param[in,out] neededVars Names of variables to try to produce
*
*/
oops::patch::Variables Vader::changeVarTraj(atlas::FieldSet & afieldset,
                      oops::patch::Variables & neededVars) {
    util::Timer timer(classname(), "changeVar");
    oops::Log::trace() << "entering Vader::changeVarTraj " << std::endl;
    oops::Log::info() << "Variables requested from Vader (Trajectory): " << neededVars.variables()
        << std::endl;

    // onlyIngredientVars are the variables that Vader should NOT later return in varsProduced
    // (i.e. variables both in afieldset AND neededVars will be returned in varsProduced)
    oops::patch::Variables onlyIngredientVars(afieldset.field_names());
    onlyIngredientVars -= neededVars;

    auto ingredientNames = afieldset.field_names();
    oops::Log::info() << "Variables provided to Vader: " << ingredientNames << std::endl;
    // Loop through all the requested fields in neededVars
    // Since neededVars can be modified by planVariable and planVariable calls
    // itself recursively, we make a copy of the list here before we start.
    std::vector<std::string> targetVariables{neededVars.variables()};
    bool recipesNeedTLAD = true;  // Only plan recipes with TL/AD methods implemented

    for (auto targetVariable : targetVariables) {
        oops::Log::debug() <<
            "Vader::changeVarTraj processing top-level requested variable: "
            << targetVariable << std::endl;
        // Since this function is called recursively, make sure targetVariable is
        // still needed
        if (!neededVars.has(targetVariable)) {
            oops::Log::debug() << targetVariable <<
                " is no longer in the variable list neededVars. Moving on." << std::endl;
            continue;
        }
        oops::patch::Variables excludedVars;
        planVariable(ingredientNames, neededVars, targetVariable, recipesNeedTLAD,
                     excludedVars, recipeExecutionPlan_);
    }
    executePlanNL(afieldset, recipeExecutionPlan_);
    // Save the trajectory in Vader's private variable
    trajectory_.clear();
    for (const auto from_Field : afieldset) {
        // Make a deep copy of each field and put it in trajectory_
        atlas::Field to_Field = from_Field.functionspace().createField<double>(
                    atlas::option::name(from_Field.name()) |
                    atlas::option::levels(from_Field.levels()));
        auto from_view = atlas::array::make_view<double, 2>(from_Field);
        auto to_view = atlas::array::make_view<double, 2>(to_Field);
        to_view.assign(from_view);
        trajectory_.add(to_Field);
    }

    if (neededVars.variables().size()) {
        oops::Log::info() << "Requested variables Vader could not produce: "
            << neededVars.variables() << std::endl;
    } else {
        oops::Log::info() << "Vader produced all requested variables."  << std::endl;
    }
    oops::patch::Variables varsProduced(afieldset.field_names());
    varsProduced -= onlyIngredientVars;
    oops::Log::trace() << "leaving Vader::changeVarTraj" << std::endl;
    return varsProduced;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Change Variable Tangent Linear
*
* \details **changeVarTL** is called externally to perform the tangent linear (TL) variable change.
* Note that unlike changeVar and changeVarTraj, the vader planVariable algorithm to determine which
* recipes to call is not invoked. Instead, the same recipe plan determined during changeVarTraj is
* executed, calling the TL methods of the planned recipes. Also note that neededVars
* is not used as input, but is kept for interface consistency, and as another way to
* communicate back to the caller which variables VADER has created.
*
* \param[in,out] afieldset This is the input/output fieldset, same as in changeVar
* \param[in,out] neededVars Names of variables to try to produce
*
*/
oops::patch::Variables Vader::changeVarTL(atlas::FieldSet & afieldset,
                      oops::patch::Variables & neededVars) const {
    oops::Log::trace() << "entering Vader::changeVarTL" << std::endl;
    oops::patch::Variables varsPopulated;
    executePlanTL(afieldset, recipeExecutionPlan_);
    for (auto varplan : recipeExecutionPlan_) {
        varsPopulated.push_back(varplan.first);
    }
    neededVars -= varsPopulated;
    oops::Log::trace() << "leaving Vader::changeVarTL" << std::endl;
    return varsPopulated;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Change Variable Adjoint
*
* \details **changeVarAD** is called externally to perform the adjoint (AD) variable change.
* Note that unlike changeVar and changeVarTraj, the vader planVariable algorithm to determine which
* recipes to call is not invoked. Instead, the same recipe plan determined during changeVarTraj is
* executed, but in reverse order, calling the AD methods of the planned recipes.
* Also note that varsToAdjoint here should be the SAME VARIABLES that are passed to changeVarTraj
* and changeVarTL in the 'neededVars' parameter. But in this case these variables should already be
* populated in afieldset. VADER will perform the adjoint methods of the recipes that produce these
* variables.
*
* \param[in,out] afieldset This is the input/output fieldset
* \param[in,out] varsToAdjoint Same vars as 'neededVars' in changeVarTraj and changeVarTL
*
*/
oops::patch::Variables Vader::changeVarAD(atlas::FieldSet & afieldset,
                      oops::patch::Variables & varsToAdjoint) const {
    oops::Log::trace() << "entering Vader::changeVarAD" << std::endl;
    oops::patch::Variables varsAdjointed;
    executePlanAD(afieldset, recipeExecutionPlan_);
    for (auto varplan : recipeExecutionPlan_) {
        oops::Log::debug() << "Adding to varsAdjointed: " << varplan.first << std::endl;
        varsAdjointed.push_back(varplan.first);
    }
    varsToAdjoint -= varsAdjointed;
    oops::Log::trace() << "leaving Vader::changeVarAD" << std::endl;
    return varsAdjointed;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Plan Variable
*
* \details **planVariable** contains Vader's primary algorithm for attempting to
* find a recipe or sequence of recipes that can produce a desired variable. It:
* * Checks the cookbook for recipes for the desired field (the targetVariable)
* * Checks each recipe to see if its required ingredients have been provided
* * If an ingredient is missing, recursively calls itself to attempt to get it
* * Adds the variable and recipe name to the "recipeExecutionPlan" if the recipe is viable.
* * If successful, removes the targetVariable from neededVars and returns 'true'
*
* \param[in] ingredientNames Names of available ingredients
* \param[in,out] neededVars Names of variables to try to produce
* \param[in] targetVariable variable name this instance is trying to plan
* \param[in] needsTLDA Flag to only consider recipes that have TLAD implemented
* \param[in,out] excludedVars variables to not try to create due to recursion chain
* \param[in,out] plan ordered list of viable recipes that will get exectued later
* \return boolean 'true' if it successfully creates a plan for targetVariable, else false
*
*/
bool Vader::planVariable(std::vector<std::string> & ingredientNames,
                         oops::patch::Variables & neededVars,
                         const std::string targetVariable,
                         const bool needsTLAD,
                         oops::patch::Variables & excludedVars,
                         vaderPlanType & plan) const {
    bool variablePlanned = false;

    oops::Log::trace() << "entering Vader::planVariable for variable: " << targetVariable <<
        std::endl;

    // Check if needed variable is part of the ingredients already.
    bool variableExists =
        (std::find(ingredientNames.begin(), ingredientNames.end(), targetVariable)
                        != ingredientNames.end());
    if (variableExists) {
        oops::Log::debug() << "Variable is part of the ingredients already. " << std::endl;
        neededVars -= targetVariable;
        return true;  // Don't need to check any recipes.
    }

    auto recipeList = cookbook_.find(targetVariable);

    // If recipeList is found, recipeList->second is a vector of unique_ptr's
    // to Recipe objects that produce 'variableName'
    if ((recipeList != cookbook_.end()) && !recipeList->second.empty()) {
        oops::Log::debug() <<
            "Vader cookbook contains at least one recipe for '" << targetVariable << "'" <<
            std::endl;
        for (const auto & recipe : recipeList->second) {
            if (needsTLAD && !recipe->hasTLAD()) {
                oops::Log::debug() << "Not checking recipe: '" << recipe->name() <<
                    "' since it does not have TL/AD methods implemented.";
                continue;
            }
            // Check that this recipe is not already in the plan to prevent
            // infinite recursion.
            for (const auto planRecipe : plan) {
                if (planRecipe.second->name() == recipe->name()) {
                    break;
                }
            }
            oops::Log::debug() << "Checking to see if we have ingredients for recipe: " <<
                recipe->name() << std::endl;
            bool haveIngredient = false;
            for (auto ingredient : recipe->ingredients()) {
                if (ingredient == targetVariable) {
                    oops::Log::error() << "Error: Ingredient list for " <<
                        recipe->name() << " contains the target." << std::endl;
                    break;
                }
                haveIngredient =
                    (std::find(ingredientNames.begin(), ingredientNames.end(), ingredient)
                        != ingredientNames.end());
                if (!haveIngredient) {
                    oops::Log::debug() << "ingredient " << ingredient << " for recipe "
                        << recipe->name() << " not in FieldSet." << std::endl;
                    if (excludedVars.has(ingredient)) {
                        oops::Log::debug() << "This ingredient is already being worked on at a " <<
                        "higher level of planVariable recursion, so not trying to create it at "
                        "this level." << std::endl;
                        haveIngredient = false;
                    } else {
                        oops::Log::debug() << "Seeing if we can create this ingredient " <<
                            "by calling planVariable recursively." << std::endl;
                        // Put this the current targetVariable in excludedVars so we don't try to
                        // create it again as an ingredient at a lower level of recursion. Then
                        // call planVariable recursively.
                        excludedVars.push_back(targetVariable);
                        haveIngredient = planVariable(ingredientNames, neededVars, ingredient,
                                                    needsTLAD, excludedVars, plan);
                        // Remove the ingredient from excludedVars since we're back to this level.
                        excludedVars -= targetVariable;
                    }
                }
                oops::Log::debug() << "ingredient " << ingredient << " for recipe "
                     << recipe->name() << (haveIngredient ? " is" : " is not") << " available."
                     << std::endl;
                if (!haveIngredient) break;  // Missing an ingredient. Don't check the others.
            }
            if (haveIngredient) {
                oops::Log::debug() <<
                    "All ingredients are in the fieldset. Adding recipe " << recipe->name() <<
                    " to recipeExecutionPlan." << std::endl;
                plan.push_back(std::pair<std::string, const std::unique_ptr<RecipeBase> &>
                                                                    ({targetVariable, recipe}));
                variablePlanned = true;
                neededVars -= targetVariable;
                ingredientNames.push_back(targetVariable);
                break;  // Found a viable recipe. Don't need to check any other potential recipes.
            } else {
                oops::Log::debug() <<
                    "Do not have and cannot create all the ingredients for recipe " <<
                    recipe->name() << std::endl;
            }
        }
    } else {
        oops::Log::debug() << "Vader cookbook does not contain a recipe for: "
            << targetVariable << std::endl;
    }
    oops::Log::trace() << "leaving Vader::planVariable for variable: " << targetVariable <<
        std::endl;
    return variablePlanned;
}

/*! \brief creates a field with correct functionspace and levels, and optionally
 *         initializes to zero. If the field already exists, asserts that it has
 *         enough levels in it.
 * \param[in,out] afieldset  Fieldset that needs to have the requested field
 *                           allocated
 * \param[in]     fieldname  Name of the field to be created or checked for dims
 * \param[in]     fs         Functionspace for creating the new field
 * \param[in]     nlevels    Number of levels for creating the new field
 * \param[in]     initZero   Initialize field with zeros? (default: false)
 */
void checkOrAddField(atlas::FieldSet & afieldset, const std::string & fieldname,
                     const atlas::FunctionSpace & fs, int nlevels,
                     bool initZero = false) {
    oops::Log::trace() << "entering checkOrAddField for Field: " << fieldname << std::endl;
    if (afieldset.has(fieldname))
    {
        // Verify the number of levels in the Field is enough for the recipe
        ASSERT(afieldset.field(fieldname).levels() >= nlevels);
    } else {
        // Create the field and put it in the FieldSet
        atlas::Field newField = fs.createField<double>(
                    atlas::option::name(fieldname) |
                    atlas::option::levels(nlevels));
        if (initZero) {
            // A new field for adjoint needs to be zeroed out
            auto view = atlas::array::make_view<double, 2>(newField);
            for (int j0 = 0; j0 < newField.shape(0); ++j0) {
                for (int j1 = 0; j1 < newField.shape(1); ++j1) {
                    view(j0, j1) = 0.0;
                }
            }
        }
        afieldset.add(newField);
    }
    oops::Log::trace() << "leaving checkOrAddField" << std::endl;
}

// ------------------------------------------------------------------------------------------------
/*! \brief Execute Plan (non-linear)
*
* \details **executePlanNL** calls, in order, the 'executeNL' method of the
* recipes specified in the recipeExecutionPlan that is passed in. (The recipeExecutionPlan is
* created through calls to planVariable.) Before executing the recipe, validations
* are performed on the passed FieldSet to ensure it has the ingredients, and the
* product Field is created and added if not already present.
*
* \param[in,out] afieldset A fieldset contaning only populated ingredient fields
* \param[in] recipeExecutionPlan ordered list of recipes that are to be exectued
*
*/
void Vader::executePlanNL(atlas::FieldSet & afieldset,
                          const vaderPlanType & recipeExecutionPlan) const {
    oops::Log::trace() << "entering Vader::executePlanNL" <<  std::endl;
    for (auto varPlan : recipeExecutionPlan) {
        for (auto ingredient :  varPlan.second->ingredients()) {
            ASSERT(afieldset.has(ingredient));
        }
        // add product field to the fieldset if needed
        checkOrAddField(afieldset, varPlan.first,
                        varPlan.second->productFunctionSpace(afieldset),
                        varPlan.second->productLevels(afieldset));
        if (varPlan.second->requiresSetup()) {
            varPlan.second->setup(afieldset);
        }
        const bool recipeSuccess = varPlan.second->executeNL(afieldset);
        ASSERT(recipeSuccess);  // At least for now, we'll require the execution to be successful
        oops::Log::info() << "Variable '" << varPlan.first <<
            "' calculated using Vader recipe " << varPlan.second->name() << std::endl;
    }
    oops::Log::trace() << "leaving Vader::executePlanNL" <<  std::endl;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Execute Plan (tangent linear)
*
* \details **executePlanTL** calls, in order, the 'execute' (tangent linear) method of the
* recipes specified in the recipeExecutionPlan that is passed in. (The recipeExecutionPlan is
* created through calls to planVariable.)
*
* \param[in,out] afieldset A fieldset contaning only populated ingredient fields
* \param[in] recipeExecutionPlan ordered list of recipes that are to be exectued
*
*/
void Vader::executePlanTL(atlas::FieldSet & afieldset,
                          const vaderPlanType & recipeExecutionPlan) const {
    oops::Log::trace() << "entering Vader::executePlanTL" <<  std::endl;
    // We must get the recipes specified in the recipeExecutionPlan out of the cookbook,
    // where they live
    for (auto varPlan : recipeExecutionPlan) {
        oops::Log::debug() << "Vader calculating variable '" << varPlan.first <<
            "' using TL recipe with name: " << varPlan.second->name() << std::endl;
        // add product field to the fieldset if needed
        checkOrAddField(afieldset, varPlan.first,
                        varPlan.second->productFunctionSpace(afieldset),
                        varPlan.second->productLevels(afieldset));
        for (auto ingredient :  varPlan.second->ingredients()) {
            ASSERT(afieldset.has(ingredient));
        }
        if (varPlan.second->requiresSetup()) {
            varPlan.second->setup(afieldset);
        }
        const bool recipeSuccess =
            varPlan.second->executeTL(afieldset, trajectory_);
        ASSERT(recipeSuccess);  // At least for now, we'll require the execution to be successful
    }
    oops::Log::trace() << "leaving Vader::executePlanTL" <<  std::endl;
}
// ------------------------------------------------------------------------------------------------
/*! \brief Execute Plan (adjoint)
*
* \details **executePlanAD** calls, in reverse order, the 'execute' (adjoint) method of the
* recipes specified in the recipeExecutionPlan that is passed in. (The recipeExecutionPlan is
* created through calls to planVariable.)
*
* \param[in,out] afieldset A fieldset contaning only populated ingredient fields
* \param[in] recipeExecutionPlan ordered list of recipes that are to be exectued
*
*/
void Vader::executePlanAD(atlas::FieldSet & afieldset,
                          const vaderPlanType & recipeExecutionPlan) const {
    oops::Log::trace() << "entering Vader::executePlanAD" <<  std::endl;
    // We must get the recipes specified in the recipeExecutionPlan out of the cookbook,
    // where they live
    // We execute the adjoints in reverse order of the recipeExecutionPlan
    for (auto varPlanIt = recipeExecutionPlan.rbegin();
         varPlanIt != recipeExecutionPlan.rend();
         ++varPlanIt) {
        oops::Log::debug()  << "Vader performing adjoint of recipe with name: " <<
            varPlanIt->second->name() << std::endl;
        ASSERT(afieldset.has(varPlanIt->first));
        // Check if ingredients need to be allocated
        for (auto ingredient :  varPlanIt->second->ingredients()) {
            // add ingredient field to the fieldset and zero out if needed
            checkOrAddField(afieldset, ingredient,
                     trajectory_.field(ingredient).functionspace(),
                     trajectory_.field(ingredient).levels(), true);
        }
        if (varPlanIt->second->requiresSetup()) {
            varPlanIt->second->setup(afieldset);
        }
        const bool recipeSuccess =
            varPlanIt->second->executeAD(afieldset, trajectory_);
        ASSERT(recipeSuccess);  // At least for now, we'll require the execution to be successful
    }
    oops::Log::trace() << "leaving Vader::executePlanAD" <<  std::endl;
}
// ------------------------------------------------------------------------------------------------
/*! \brief print
*
* \details **print** this method prints the current Vader cookbook in yaml form to an osstream
*
* \param[in] os the stream to print to
*
*/
void Vader::print(std::ostream & os) const {
    oops::Log::trace() << "entering Vader::print" <<  std::endl;
    std::string cookbookString("VADER Cookbook\ncookbook:\n");
    for (const auto & cbMember : cookbook_) {
        cookbookString += "  ";  // 2 space indent for yaml
        bool firstRecipe = true;
        cookbookString += cbMember.first + ": [";
        for (const auto & recipe : cbMember.second) {
            if (firstRecipe) {
                firstRecipe = false;
            } else {
                cookbookString += ",";
            }
            cookbookString += recipe->name();
        }
        cookbookString += "]\n";
    }
    os << cookbookString;
    oops::Log::trace() << "leaving Vader::print" <<  std::endl;
}
}  // namespace vader
