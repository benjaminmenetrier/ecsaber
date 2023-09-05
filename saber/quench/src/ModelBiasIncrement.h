/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>

#include "util/Printable.h"

#include "src/Geometry.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class ModelBias;
  class ModelBiasCovariance;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasIncrement : public util::Printable {
 public:
/// Constructor, destructor
  ModelBiasIncrement(const Geometry &, const eckit::Configuration &);
  ModelBiasIncrement(const ModelBiasIncrement &, const bool = true);
  ModelBiasIncrement(const ModelBiasIncrement &, const eckit::Configuration &);
  ~ModelBiasIncrement();

/// Linear algebra operators
  void diff(const ModelBias &, const ModelBias &) {}
  void zero() {}
  ModelBiasIncrement & operator=(const ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator+=(const ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator-=(const ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator*=(const double &) {return *this;}
  void axpy(const double, const ModelBiasIncrement &) {}
  double dot_product_with(const ModelBiasIncrement &) const {return 0.0;}

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
