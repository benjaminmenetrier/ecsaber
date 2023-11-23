/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SqrtOfSpectralCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SqrtOfSpectralCovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbReadParameters> readParams{"read", this};
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class SqrtOfSpectralCovariance : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SqrtOfSpectralCovariance";}

  typedef SqrtOfSpectralCovarianceParameters Parameters_;

  SqrtOfSpectralCovariance(const oops::GeometryData &,
                           const oops::patch::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);

  virtual ~SqrtOfSpectralCovariance() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return outerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;

  void read() override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::patch::Variables activeVars_;
  /// Outer variables
  oops::patch::Variables outerVars_;

  /// Covariance statistics
  atlas::FieldSet spectralUMatrices_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Geometry data
  const oops::GeometryData & innerGeometryData_;
};

}  // namespace spectralb
}  // namespace saber
