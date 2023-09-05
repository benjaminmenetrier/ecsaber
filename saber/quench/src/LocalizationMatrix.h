/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"

#include "util/abor1_cpp.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace quench {
  class Geometry;
  class Increment;
  class IncrEnsCtlVec;

/// Localization matrix for quench model.

// -----------------------------------------------------------------------------
class LocalizationMatrix: public util::Printable,
                            private eckit::NonCopyable,
                            private util::ObjectCounter<LocalizationMatrix> {
 public:
  static const std::string classname() {return "quench::LocalizationMatrix";}

  LocalizationMatrix(const Geometry &, const eckit::Configuration &)
    {ABORT("not implemented yet");}
  ~LocalizationMatrix() {}
  void multiply(Increment &) const
    {ABORT("not implemented yet");}
  void multiplySqrt(const IncrEnsCtlVec &, Increment &) const
    {ABORT("not implemented yet");}
  void multiplySqrtTrans(const Increment &, IncrEnsCtlVec &) const
    {ABORT("not implemented yet");}

 private:
  void print(std::ostream &) const
    {ABORT("not implemented yet");}
};
// -----------------------------------------------------------------------------

}  // namespace quench
