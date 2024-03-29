/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace oops {
  namespace patch{
class Variables;
}
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class HydroBalParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydroBalParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::patch::Variables mandatoryActiveVars() const override {
    return oops::patch::Variables({
        "air_pressure_levels",
        "hydrostatic_exner_levels",
        "virtual_potential_temperature"});
  }

  oops::patch::Variables activeInnerVars(const oops::patch::Variables& outerVars) const override {
    oops::patch::Variables vars({"air_pressure_levels",
                          "hydrostatic_exner_levels"});
    const int modelLevels = outerVars.getLevels("virtual_potential_temperature");
    vars.addMetaData("air_pressure_levels", "levels", modelLevels + 1);
    vars.addMetaData("hydrostatic_exner_levels", "levels", modelLevels + 1);
    return vars;
  }

  oops::patch::Variables activeOuterVars(const oops::patch::Variables& outerVars) const override {
    oops::patch::Variables vars({"virtual_potential_temperature"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

class HydroBal : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydroBal";}

  typedef HydroBalParameters Parameters_;

  HydroBal(const oops::GeometryData &,
           const oops::patch::Variables &,
           const eckit::Configuration &,
           const Parameters_ &,
           const oops::FieldSet3D &,
           const oops::FieldSet3D &);
  virtual ~HydroBal();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::patch::Variables innerVars_;
  const oops::patch::Variables activeOuterVars_;
  const oops::patch::Variables innerOnlyVars_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
