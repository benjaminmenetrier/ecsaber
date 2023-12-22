/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <string>

#include "oops/interface/GenericMatrix.h"

#include "src/Traits.h"

#include "util/abor1_cpp.h"
#include "util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Increment;
  class State;

/// Hybrid weight class: Apply hybrid weight

// -----------------------------------------------------------------------------

class HybridWeight : public oops::GenericMatrix<Traits>,
                     private util::ObjectCounter<HybridWeight> {
 public:
  static const std::string classname() {return "quench::HybridWeight";}

  explicit HybridWeight(const eckit::Configuration &);
  ~HybridWeight();

  void multiply(State &) const
    {ABORT("not implemented yet");}
  void multiply(Increment &) const;
  void inverseMultiply(Increment &) const
    {ABORT("not implemented yet");}
  void transposeInverseMultiply(Increment &) const
    {ABORT("not implemented yet");}

 private:
  void print(std::ostream &) const;

  double wgt_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
