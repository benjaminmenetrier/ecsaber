/*
 * (C) Copyright 2023 UCAR
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
#include "OceanConversions/OceanConversions.interface.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader {

// -------------------------------------------------------------------------------------------------

class SeaWaterTemperature_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(SeaWaterTemperature_AParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{"recipe name", this};
};

/*! \brief SeaWaterTemperature_A class defines a recipe for insitu water temperature
 *
 *  \details This instantiation of RecipeBase produces insitu water temperature using
 *           sea water potential temperature, salinity, longitude, lattitude, and
 *           sea_water_depth information.
 */
class SeaWaterTemperature_A : public RecipeBase {
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef SeaWaterTemperature_AParameters Parameters_;

    SeaWaterTemperature_A(const Parameters_ &, const VaderConfigVars &);

    // Recipe base class overrides
    std::string name() const override;
    std::string product() const override;
    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
    const VaderConfigVars & configVariables_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace vader
