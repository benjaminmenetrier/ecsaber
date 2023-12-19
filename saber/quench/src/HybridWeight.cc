/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/HybridWeight.h"

#include <cmath>
#include <vector>

#include "util/Logger.h"

#include "src/Increment.h"

namespace quench {

// -----------------------------------------------------------------------------
HybridWeight::HybridWeight(const eckit::Configuration & config) :
  wgt_(std::sqrt(config.getDouble("weight")))
{
  oops::Log::trace() << "HybridWeight constructed" << std::endl;
}
// -----------------------------------------------------------------------------
HybridWeight::~HybridWeight() {
  oops::Log::trace() << "HybridWeight destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void HybridWeight::multiply(Increment & dx) const {
  dx *= wgt_;
}
// -----------------------------------------------------------------------------
void HybridWeight::print(std::ostream & os) const {
  os << std::endl << "Hybrid weight:" << wgt_;
}
// -----------------------------------------------------------------------------

}  // namespace quench
