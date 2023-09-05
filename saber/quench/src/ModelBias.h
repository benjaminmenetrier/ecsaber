/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "eckit/memory/NonCopyable.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Model;
  class ModelBiasIncrement;

/// Model error for the  model.
/*!
 * This class is used to manipulate parameters of the model that
 * can be estimated in the assimilation. This includes model bias for
 * example but could be used for other parameters to be estimated.
 * This is sometimes referred to as augmented state or augmented
 * control variable in the litterature.
 * The augmented state is understood here as an augmented 4D state.
 */

// -----------------------------------------------------------------------------

class ModelBias : public util::Printable,
                  private eckit::NonCopyable,
                  private util::ObjectCounter<ModelBias> {
 public:
  static const std::string classname() {return "quench::ModelBias";}

  ModelBias(const Geometry &, const Model &,
            const eckit::Configuration &);
  ModelBias(const Geometry &, const eckit::Configuration &);
  ModelBias(const Geometry &, const ModelBias &);
  ModelBias(const ModelBias &, const bool);
  ~ModelBias();

  ModelBias & operator+=(const ModelBiasIncrement &) {return *this;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

/// Status
  const bool & isActive() const {return active_;}

 private:
  void print(std::ostream & os) const {}
  const bool active_ = false;
};

// -----------------------------------------------------------------------------

}  // namespace quench
