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
#include "saber/spectralb/CovarianceStatistics.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SqrtOfMOSpectralCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SqrtOfMOSpectralCovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbParameters> readParams{"read", this};
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class SqrtOfMOSpectralCovariance : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SqrtOfMOSpectralCovariance";}

  typedef SqrtOfMOSpectralCovarianceParameters Parameters_;

  SqrtOfMOSpectralCovariance(const oops::GeometryData &,
                           const oops::patch::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);

  virtual ~SqrtOfMOSpectralCovariance() = default;

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

  /// Option to use sqrt of vertical covariances or correlations
  bool variance_opt_;
  /// Covariance statistics
  // Note: only need square root of vertical covariances or correlations from this;
  // probably can be gotten in the ctor and saved here instead of cs_
  std::unique_ptr<CovStat_ErrorCov> cs_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Geometry data
  const oops::GeometryData & innerGeometryData_;
};

}  // namespace spectralb
}  // namespace saber
