/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "util/DateTime.h"
#include "util/dot_product.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;  // TODO(Benjamin): for cy46 only, should be removed later
  class Covariance;

  /// IncrModCtlVec Class
  /*!
   *  Increment part of the control vector (acted upon by the static background
   *  error covariance matrix)
   */

// -----------------------------------------------------------------------------

class IncrModCtlVec : public util::Printable,
                      private util::ObjectCounter<IncrModCtlVec> {
 public:
  static const std::string classname() {return "quench::IncrModCtlVec";}

/// Constructor, destructor
  IncrModCtlVec();
  explicit IncrModCtlVec(const Covariance &);
  IncrModCtlVec(const Covariance &, const IncrModCtlVec &);
  IncrModCtlVec(const IncrModCtlVec &, const bool = true);
  virtual ~IncrModCtlVec();

  // TODO(Benjamin): for cy46 only, should be removed later
  IncrModCtlVec(const Geometry &, const Covariance &);
  IncrModCtlVec(const Geometry &, const Covariance &, const IncrModCtlVec &);

/// Basic operators
  void zero() {}
  IncrModCtlVec & operator =(const IncrModCtlVec &) {return *this;}
  IncrModCtlVec & operator+=(const IncrModCtlVec &) {return *this;}
  IncrModCtlVec & operator-=(const IncrModCtlVec &) {return *this;}
  IncrModCtlVec & operator*=(const double &) {return *this;}
  void axpy(const double &, const IncrModCtlVec &, const bool check = true) {}
  double dot_product_with(const IncrModCtlVec &) const {return 0.0;}
  void schur_product_with(const IncrModCtlVec &) {}
  void random() {}

/// Section
  void getSection(std::vector<double> &,
                 const int &,
                 const int &) const {}
  void setSection(const std::vector<double> &,
                 const int &,
                 const int &) {}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

/// Data
 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
