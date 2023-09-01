/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <string>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Model;
  class ModelBias;
  class State;

/// Model error neural network model for the  model.
/*!
 * This class is ..
 */

// -----------------------------------------------------------------------------

class ModelBiasEstimator : public util::Printable,
                             private eckit::NonCopyable,
                             private util::ObjectCounter<ModelBiasEstimator> {
 public:
  static const std::string classname() {return "quench::ModelBiasEstimator";}

/// Constructors, destructors
  ModelBiasEstimator(const Geometry &,
                     const Model &,
                     const eckit::Configuration &) {}
  ~ModelBiasEstimator() {}

/// Methods
  void estimate(const State &, ModelBias &) const {}
  void estimate(const State &, const State &, ModelBias &) const {}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace quench
