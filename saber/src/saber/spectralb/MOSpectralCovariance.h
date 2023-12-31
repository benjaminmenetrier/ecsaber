/*
 * (C) Crown Copyright 2022-2023 Met Office
 * (C) Copyright 2022- UCAR
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
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/spectralb/CovarianceStatistics.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class MOSpectralCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MOSpectralCovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbParameters> readParams{"read", this};
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class MOSpectralCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::MOSpectralCovariance";}

  typedef MOSpectralCovarianceParameters Parameters_;

  MOSpectralCovariance(const oops::GeometryData &,
                     const oops::patch::Variables &,
                     const eckit::Configuration &,
                     const Parameters_ &,
                     const oops::FieldSet3D &,
                     const oops::FieldSet3D &);

  virtual ~MOSpectralCovariance() = default;

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  void read() override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::patch::Variables activeVars_;
  /// Option to use vertical covariances or correlations
  bool variance_opt_;
  /// Covariance statistics
  // Note: only need vertical covariances or correlations from this;
  // probably can be gotten in the ctor and saved here instead of cs_
  std::unique_ptr<CovStat_ErrorCov> cs_;
  /// Geometry data
  const oops::GeometryData & geometryData_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
};

}  // namespace spectralb
}  // namespace saber
