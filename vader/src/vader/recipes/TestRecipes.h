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
#include "atlas/functionspace/PointCloud.h"
#include "atlas/util/Point.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

class TestRecipeParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(TestRecipeParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{
     "recipe name",
     this};
};

// ------------------------------------------------------------------------------------------------
/*! \brief Test_VarA_from_B class defines a recipe for TestVarA
 *
 *  \details This instantiation of RecipeBase produces TestVarA
 *           using TestVarB as input. Non-linear recipe only.
 */
class Test_VarA_from_B : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef TestRecipeParameters Parameters_;

    Test_VarA_from_B(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override {return "TestVarA";}
    std::vector<std::string> ingredients() const override {return {"TestVarB"};}
    size_t productLevels(const atlas::FieldSet &) const override {return 0;}
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override
                                {return atlas::functionspace::PointCloud({
                                    atlas::PointXY(0.0, 0.0)});}
    // bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    // bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    // bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};

// ------------------------------------------------------------------------------------------------
/*! \brief Test_VarA_from_C_D class defines a recipe for TestVarA
 *
 *  \details This instantiation of RecipeBase produces TestVarA
 *           using TestVarC and TestVarD as input. Non-linear recipe only.
 */
class Test_VarA_from_C_D : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef TestRecipeParameters Parameters_;

    Test_VarA_from_C_D(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override {return "TestVarA";}
    std::vector<std::string> ingredients() const override {return {"TestVarC", "TestVarD"};}
    size_t productLevels(const atlas::FieldSet &) const override {return 0;}
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override
                                {return atlas::functionspace::PointCloud({
                                    atlas::PointXY(0.0, 0.0)});}
    // bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    // bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    // bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};

// ------------------------------------------------------------------------------------------------
/*! \brief Test_VarB_from_E class defines a recipe for TestVarB
 *
 *  \details This instantiation of RecipeBase produces TestVarB
 *           using TestVarE as input. Non-linear recipe only.
 */
class Test_VarB_from_E : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef TestRecipeParameters Parameters_;

    Test_VarB_from_E(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override {return "TestVarB";}
    std::vector<std::string> ingredients() const override {return {"TestVarE"};}
    size_t productLevels(const atlas::FieldSet &) const override {return 0;}
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override
                                {return atlas::functionspace::PointCloud({
                                    atlas::PointXY(0.0, 0.0)});}
    // bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    // bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    // bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};
// ------------------------------------------------------------------------------------------------
/*! \brief Test_VarB_from_A class defines a recipe for TestVarB
 *
 *  \details This instantiation of RecipeBase produces TestVarB
 *           using TestVarA as input. Non-linear recipe only.
 *           This recipe combines with recipe Test_VarA_from_B to make 
 *           sure the Vader algorithm does not end up in an infinte
 *           recursion loop when attempting to create ingredients.
 */
class Test_VarB_from_A : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef TestRecipeParameters Parameters_;

    Test_VarB_from_A(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override {return "TestVarB";}
    std::vector<std::string> ingredients() const override {return {"TestVarA"};}
    size_t productLevels(const atlas::FieldSet &) const override {return 0;}
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override
                                {return atlas::functionspace::PointCloud({
                                    atlas::PointXY(0.0, 0.0)});}
    // bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    // bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    // bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};
// ------------------------------------------------------------------------------------------------
/*! \brief Test_VarC_from_F class defines a recipe for TestVarC
 *
 *  \details This instantiation of RecipeBase produces TestVarC
 *           using TestVarF as input. Non-linear recipe only.
 */
class Test_VarC_from_F : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef TestRecipeParameters Parameters_;

    Test_VarC_from_F(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override {return "TestVarC";}
    std::vector<std::string> ingredients() const override {return {"TestVarF"};}
    size_t productLevels(const atlas::FieldSet &) const override {return 0;}
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override
                                {return atlas::functionspace::PointCloud({
                                    atlas::PointXY(0.0, 0.0)});}
    // bool hasTLAD() const override { return true; }
    bool executeNL(atlas::FieldSet &) override;
    // bool executeTL(atlas::FieldSet &, const atlas::FieldSet &) override;
    // bool executeAD(atlas::FieldSet &, const atlas::FieldSet &) override;

 private:
};

}  // namespace vader
