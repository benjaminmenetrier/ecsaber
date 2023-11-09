/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Geometry.h"

// Forward declarations
namespace quench {
  class Increment;
  class IncrModCtlVec;
  class State;
  class Variables;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for quench model.

class Covariance : public util::Printable,
                   private eckit::NonCopyable,
                   private util::ObjectCounter<Covariance> {
 public:
  static const std::string classname() {return "quench::Covariance";}

  Covariance(const Geometry &, const Variables &,
             const eckit::Configuration &, const State &)
    {throw eckit::NotImplemented(Here());}
  ~Covariance() {}

  void linearize(const State &, const Geometry &)
    {throw eckit::NotImplemented(Here());}
  void multiply(const Increment &, Increment &) const
    {throw eckit::NotImplemented(Here());}
  void inverseMultiply(const Increment &, Increment &) const
    {throw eckit::NotImplemented(Here());}
  void multiplySqrt(const IncrModCtlVec &, Increment &) const
    {throw eckit::NotImplemented(Here());}
  void multiplySqrtTrans(const Increment &, IncrModCtlVec &) const
    {throw eckit::NotImplemented(Here());}
  void randomize(Increment &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream & os) const {os << "Covariance";}
};
// -----------------------------------------------------------------------------

}  // namespace quench
