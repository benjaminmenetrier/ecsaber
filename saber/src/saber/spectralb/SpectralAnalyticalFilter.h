/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralAnalyticalFilterParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralAnalyticalFilterParameters, SaberBlockParametersBase)

 public:
    /// Whether to normalize as a localization function
    oops::Parameter<bool> normalize{"normalize filter variance", true, this};
    /// Filter specifications
    oops::Parameter<eckit::LocalConfiguration> function{"function",
                                                        eckit::LocalConfiguration(), this};

    oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class SpectralAnalyticalFilter : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralAnalyticalFilter";}

  typedef SpectralAnalyticalFilterParameters Parameters_;

  SpectralAnalyticalFilter(const oops::GeometryData &,
                           const oops::patch::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);

  virtual ~SpectralAnalyticalFilter() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                         const oops::patch::Variables & innerVars) const override;

  oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                         const oops::patch::Variables & outerVars) const override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::patch::Variables activeVars_;
  /// inner Geometry Data for next block
  const oops::GeometryData & innerGeometryData_;
  /// inner variables for next block
  const oops::patch::Variables innerVars_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Filter in spectral space
  const std::vector<double> spectralFilter_;
};

}  // namespace spectralb
}  // namespace saber
