/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"

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
  IncrEnsCtlVec()
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  IncrEnsCtlVec(const Geometry &, const LocalizationMatrix &)
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  IncrEnsCtlVec(const Geometry &, const LocalizationMatrix &, const IncrEnsCtlVec &)
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  IncrEnsCtlVec(const IncrEnsCtlVec &, const bool)
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  virtual ~IncrEnsCtlVec() {}

/// Basic operators
  void zero()
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  IncrEnsCtlVec & operator =(const IncrEnsCtlVec &) {return *this;}
  IncrEnsCtlVec & operator+=(const IncrEnsCtlVec &) {return *this;}
  IncrEnsCtlVec & operator-=(const IncrEnsCtlVec &) {return *this;}
  IncrEnsCtlVec & operator*=(const double &) {return *this;}
  void axpy(const double &, const IncrEnsCtlVec &, const bool check = true)
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  double dot_product_with(const IncrEnsCtlVec &) const
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here()); return 0.0;}
  void schur_product_with(const IncrEnsCtlVec &)
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  void random()
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}

/// I/O and diagnostics
  void read(const eckit::Configuration &)
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  void write(const eckit::Configuration &) const
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
  double norm() const
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here()); return 0.0;}

 private:
  void print(std::ostream &) const
    {throw eckit::Exception("no IncrEnsCtlVec with quench", Here());}
};
// -----------------------------------------------------------------------------

}  // namespace quench
