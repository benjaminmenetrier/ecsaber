/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <map>
#include <ostream>
#include <string>

#include "util/DateTime.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;
  class Gom;
  class ObsVec;

// -----------------------------------------------------------------------------
/// ObsSpace empty class

class ObsSpace : public util::Printable {
 public:
  ObsSpace(const eckit::Configuration &,
             const Geometry &,
             const util::DateTime &,
             const util::DateTime &,
             const bool lscreend = false) {}
  ~ObsSpace() {}

  void screenObservations(const ObsVec &, const Gom &) const {}

  void saveObservations() const {}

  void generateDistribution(const eckit::Configuration & conf) {}

  void printJo(const ObsVec &, const ObsVec &) {}

 private:
  void print(std::ostream &) const {}
};

}  // namespace quench

