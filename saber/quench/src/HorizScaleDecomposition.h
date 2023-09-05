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

#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace quench {
  class Increment;

/// HorizScaleDecomposition matrix for quench model.

// -----------------------------------------------------------------------------
class HorizScaleDecomposition: public util::Printable,
                               private eckit::NonCopyable,
                               private util::ObjectCounter<HorizScaleDecomposition> {
 public:
    explicit HorizScaleDecomposition(const eckit::Configuration &) {}

    void collect(const Increment &, const int &) const {}
    void decompose(const Increment &, const int &) const {}

    static const std::string classname() {return "quench::HorizScaleDecomposition";}

 private:
    void print(std::ostream & os) const {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
