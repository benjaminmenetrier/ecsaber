/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "util/DateTime.h"
#include "util/dot_product.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Geometry.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class LocalizationMatrix;

  /// IncrEnsCtlVec Class
  /*!
   *  Augmented part of the control vector
   */

// -----------------------------------------------------------------------------

class IncrEnsCtlVec : public util::Printable,
                      private util::ObjectCounter<IncrEnsCtlVec> {
 public:
  static const std::string classname() {return "quench::IncrEnsCtlVec";}

/// Constructor, destructor
  IncrEnsCtlVec() {}
  IncrEnsCtlVec(const Geometry &, const LocalizationMatrix &) {}
  IncrEnsCtlVec(const Geometry &, const LocalizationMatrix &, const IncrEnsCtlVec &) {}
  IncrEnsCtlVec(const IncrEnsCtlVec &, const bool) {}
  virtual ~IncrEnsCtlVec() {}

/// Basic operators
  void zero() {}
  // void zero(const util::DateTime &);
  IncrEnsCtlVec & operator =(const IncrEnsCtlVec &) {return *this;}
  IncrEnsCtlVec & operator+=(const IncrEnsCtlVec &) {return *this;}
  IncrEnsCtlVec & operator-=(const IncrEnsCtlVec &) {return *this;}
  IncrEnsCtlVec & operator*=(const double &) {return *this;}
  void axpy(const double &, const IncrEnsCtlVec &, const bool check = true) {}
  double dot_product_with(const IncrEnsCtlVec &) const {return 0.0;}
  void schur_product_with(const IncrEnsCtlVec &) {}
  void random();

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
