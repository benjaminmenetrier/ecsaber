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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace vader
{

class SurfaceAirPressure_AParameters : public RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(SurfaceAirPressure_AParameters, RecipeParametersBase)

 public:
    oops::RequiredParameter<std::string> name{"recipe name", this};
};

// -------------------------------------------------------------------------------------------------
/// Recipe base class

/*! \brief SurfaceAirPressure_A class defines a recipe for surface pressure
 *
 *  \details This recipe produces surface pressure from air pressure thickness (delp) by summing
 *           the pressure at the model top with all the delp values. It does not provide
 *           TL/AD algorithms.
 */
class SurfaceAirPressure_A : public RecipeBase
{
 public:
    static const char Name[];
    static const std::vector<std::string> Ingredients;

    typedef SurfaceAirPressure_AParameters Parameters_;

    SurfaceAirPressure_A(const Parameters_ &, const VaderConfigVars &);

    std::string name() const override;
    std::string product() const override;
    std::vector<std::string> ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override;
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    bool executeNL(atlas::FieldSet &) override;

 private:
    std::map<std::string, double> p0Defaults_;
    const VaderConfigVars & configVariables_;
};
}  // namespace vader
