/*
 * (C) Crown Copyright 2023 Met Office
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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

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

// -----------------------------------------------------------------------------

class HpHexnerToPExnerParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HpHexnerToPExnerParameters, SaberBlockParametersBase)
 public:
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables({
    "air_pressure_levels",
    "exner_levels_minus_one",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels"});}
};

// -----------------------------------------------------------------------------
/// \brief This saber block is here to copy hydrostatic_exner_levels into
///        exner_levels_minus_one and hydrostatic_pressure_levels into
///        air_pressure_levels

class HpHexnerToPExner : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HpHexnerToPExner";}

  typedef HpHexnerToPExnerParameters Parameters_;

  HpHexnerToPExner(const oops::GeometryData &,
                   const oops::patch::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~HpHexnerToPExner();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::patch::Variables innerVars_;
  oops::patch::Variables activeVars_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
