/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>

#include "oops/interface/LinearModelBase.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace quench {
// -----------------------------------------------------------------------------
///  linear model definition.
/*!
 *   linear model definition and configuration parameters.
 */

class Tlm: public oops::LinearModelBase<Traits>,
           private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname() {return "quench::Tlm";}

  Tlm(const Geometry &, const eckit::Configuration &) {}
  ~Tlm() {}

/// Model trajectory computation
  void setTrajectory(State &, const Model &, const ModelBias &) override {}

/// Run TLM and its adjoint
  void initializeTL(Increment &, const ModelBiasIncrement &) const override {}
  void stepTL(Increment &, const ModelBiasIncrement &) const override {}
  void finalizeTL(Increment &) const override {}

  void initializeAD(Increment &, const ModelBiasIncrement &) const override {}
  void stepAD(Increment &, ModelBiasIncrement &) const override {}
  void finalizeAD(Increment &) const override {}

/// Other utilities
  const util::Duration & timeResolution() const override {}

 private:
  void print(std::ostream &) const override {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
