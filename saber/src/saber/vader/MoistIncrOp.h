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

#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/MoistIncrOpParameters.h"

namespace oops {
  namespace patch{
class Variables;
}
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class MoistIncrOp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::MoistIncrOp";}

  typedef MoistIncrOpParameters Parameters_;

  MoistIncrOp(const oops::GeometryData &,
              const oops::patch::Variables &,
              const eckit::Configuration &,
              const Parameters_ &,
              const oops::FieldSet3D &,
              const oops::FieldSet3D &);
  virtual ~MoistIncrOp();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::patch::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet & fset) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::patch::Variables innerVars_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
