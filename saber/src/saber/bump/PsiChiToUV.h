/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/bump/BUMPParameters.h"
#include "saber/bump/lib/BUMP.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class PsiChiToUVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(PsiChiToUVParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<BUMPParameters> readParams{"read", this};
  oops::OptionalParameter<BUMPParameters> calibrationParams{"calibration", this};

  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables({
    "stream_function",
    "velocity_potential",
    "eastward_wind",
    "northward_wind"});}
};

// -----------------------------------------------------------------------------

class PsiChiToUV : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::bump::PsiChiToUV";}

  typedef PsiChiToUVParameters Parameters_;

  PsiChiToUV(const oops::GeometryData &,
             const oops::patch::Variables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const oops::FieldSet3D &,
             const oops::FieldSet3D &);

  virtual ~PsiChiToUV();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

  void read() override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::patch::Variables innerVars_;
  oops::patch::Variables outerVars_;
  BUMPParameters bumpParams_;
  std::unique_ptr<bump_lib::BUMP> bump_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
