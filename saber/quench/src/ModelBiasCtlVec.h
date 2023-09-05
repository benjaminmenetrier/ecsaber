/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class ModelBias;
  class ModelBiasCovariance;

  /// ModelBiasCtlVec Class
  /*!
   *  Part of the control vector corresponding to the ModelAuxIncrement part of the
   *  control increment
   */

// -----------------------------------------------------------------------------

class ModelBiasCtlVec : public util::Printable {
 public:
/// Constructor, destructor
  explicit ModelBiasCtlVec(const ModelBiasCovariance &) {}
  ModelBiasCtlVec(const ModelBiasCovariance &, const ModelBiasCtlVec &) {}
  ModelBiasCtlVec(const ModelBiasCtlVec &, const bool copy = true) {}
  ~ModelBiasCtlVec() {}

/// Linear algebra operators
  void zero() {}
  ModelBiasCtlVec & operator=(const ModelBiasCtlVec &) {return *this;}
  ModelBiasCtlVec & operator+=(const ModelBiasCtlVec &) {return *this;}
  ModelBiasCtlVec & operator-=(const ModelBiasCtlVec &) {return *this;}
  ModelBiasCtlVec & operator*=(const double &) {return *this;}
  void axpy(const double &, const ModelBiasCtlVec &) {}
  double dot_product_with(const ModelBiasCtlVec &) const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

/// Status
  const bool & isActive() const {return active_;}

 private:
  void print(std::ostream &) const;
  const bool active_ = false;
};

// -----------------------------------------------------------------------------

}  // namespace quench
