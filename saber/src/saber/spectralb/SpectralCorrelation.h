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
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralCorrelationParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralCorrelationParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbCalibrationVertCovParameters>
    calibrationParams{"calibration", this};
  oops::OptionalParameter<spectralbReadParameters> readParams{"read", this};
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class SpectralCorrelation : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralCorrelation";}

  typedef SpectralCorrelationParameters Parameters_;

  SpectralCorrelation(const oops::GeometryData &,
                      const oops::patch::Variables &,
                      const eckit::Configuration &,
                      const Parameters_ &,
                      const oops::FieldSet3D &,
                      const oops::FieldSet3D &);

  virtual ~SpectralCorrelation() = default;

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void write() const override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::patch::Variables activeVars_;
  /// Vertical Spectral Correlations
  atlas::FieldSet spectralVerticalCorrelations_;
  /// Geometry data
  const oops::GeometryData & geometryData_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
};

}  // namespace spectralb
}  // namespace saber
