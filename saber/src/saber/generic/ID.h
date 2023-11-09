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

#include "oops/base/GeometryData.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class IDParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(IDParameters, SaberBlockParametersBase)
 public:
  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class IDCentral : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::IDCentral";}

  typedef IDParameters Parameters_;

  IDCentral(const oops::GeometryData &,
            const oops::patch::Variables &,
            const eckit::Configuration &,
            const Parameters_ &,
            const oops::FieldSet3D &,
            const oops::FieldSet3D &);

  virtual ~IDCentral();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  size_t ctlVecSize() const override {return ctlVecSize_;}
  void multiplySqrt(const atlas::Field &, atlas::FieldSet &, const size_t &) const override;
  void multiplySqrtAD(const atlas::FieldSet &, atlas::Field &, const size_t &) const override;

 private:
  const oops::GeometryData & geometryData_;
  const oops::patch::Variables activeVars_;
  size_t ctlVecSize_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

class IDOuter : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::IDOuter";}

  typedef IDParameters Parameters_;

  IDOuter(const oops::GeometryData &,
          const oops::patch::Variables &,
          const eckit::Configuration &,
          const Parameters_ &,
          const oops::FieldSet3D &,
          const oops::FieldSet3D &);

  virtual ~IDOuter() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::patch::Variables innerVars_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
