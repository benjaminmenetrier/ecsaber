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

class DryAirDensityParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryAirDensityParameters, SaberBlockParametersBase)
 public:
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables({
    "dry_air_density_levels_minus_one",
    "air_pressure_levels",
    "virtual_potential_temperature"});}

  oops::patch::Variables activeInnerVars(const oops::patch::Variables& outerVars) const override {
    oops::patch::Variables vars({"air_pressure_levels",
                          "virtual_potential_temperature"});
    const int modelLevels = outerVars.getLevels("dry_air_density_levels_minus_one");
    vars.addMetaData("air_pressure_levels", "levels", modelLevels + 1);
    vars.addMetaData("virtual_potential_temperature", "levels", modelLevels);
    return vars;
  }

  oops::patch::Variables activeOuterVars(const oops::patch::Variables& outerVars) const override {
    oops::patch::Variables vars({"dry_air_density_levels_minus_one"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

class DryAirDensity : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::DryAirDensity";}

  typedef DryAirDensityParameters Parameters_;

  DryAirDensity(const oops::GeometryData &,
                const oops::patch::Variables &,
                const eckit::Configuration &,
                const Parameters_ &,
                const oops::FieldSet3D &,
                const oops::FieldSet3D &);
  virtual ~DryAirDensity();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::patch::Variables innerVars_;
  const oops::patch::Variables activeOuterVars_;
  const oops::patch::Variables innerOnlyVars_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
