/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"

#include "util/abor1_cpp.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Geometry.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace quench {
  class Increment;
  class IncrModCtlVec;
  class State;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for quench model.

class ErrorCovariance : public util::Printable,
                        private eckit::NonCopyable,
                        private util::ObjectCounter<ErrorCovariance> {
 public:
  static const std::string classname() {return "quench::ErrorCovariance";}

  ErrorCovariance(const Geometry &, const oops::Variables &,
                  const eckit::Configuration &, const State &)
    {ABORT("not implemented yet");}
  ~ErrorCovariance() {}

  void linearize(const State &, const Geometry &)
    {ABORT("not implemented yet");}
  void multiply(const Increment &, Increment &) const
    {ABORT("not implemented yet");}
  void inverseMultiply(const Increment &, Increment &) const
    {ABORT("not implemented yet");}
  void multiplySqrt(const IncrModCtlVec &, Increment &) const
    {ABORT("not implemented yet");}
  void multiplySqrtTrans(const Increment &, IncrModCtlVec &) const
    {ABORT("not implemented yet");}
  void randomize(Increment &) const
    {ABORT("not implemented yet");}

 private:
  void print(std::ostream & os) const {os << "ErrorCovariance";}
};
// -----------------------------------------------------------------------------

}  // namespace quench
