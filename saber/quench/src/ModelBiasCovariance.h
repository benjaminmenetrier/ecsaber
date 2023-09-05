/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/memory/NonCopyable.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class ModelBias;
  class ModelBiasCtlVec;
  class ModelBiasIncrement;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasCovariance : public util::Printable,
                            private eckit::NonCopyable,
                            private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname() {return "quench::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration &, const Geometry &) {}
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const Geometry &) {}
  void multiply(const ModelBiasIncrement &, ModelBiasIncrement &) const {}
  void inverseMultiply(const ModelBiasIncrement &, ModelBiasIncrement &) const {}
  void multiplySqrt(const ModelBiasCtlVec &, ModelBiasIncrement &) const {}
  void multiplySqrtTrans(const ModelBiasIncrement &, ModelBiasCtlVec &) const {}
  void randomize(ModelBiasIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

// Status
  const bool & isActive() const {return active_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_ = eckit::LocalConfiguration();
  const bool active_ = false;
};

// -----------------------------------------------------------------------------

}  // namespace quench
