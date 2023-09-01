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
namespace oops {
  class Variables;
}

namespace quench {
  class Geometry;
  class State;
  class Increment;

/// Interpolator matrix for quench model.

// -----------------------------------------------------------------------------
class Interpolator: public util::Printable,
                            private eckit::NonCopyable,
                            private util::ObjectCounter<Interpolator> {
 public:
    Interpolator(const State &, const eckit::Configuration &, const bool) {}
    Interpolator(const State &, const Geometry &, const eckit::Configuration &, const bool) {}
    explicit Interpolator(const Interpolator &) {}
    ~Interpolator() {}

    Interpolator & operator=(const Interpolator &) {return *this;}
    Interpolator & operator+=(const State &) {return *this;}
    Interpolator & operator-=(const State &) {return *this;}

    void set(const Geometry &, const Geometry &) {}
    virtual void set_vector(const eckit::Configuration &) {}

    void interpolate(const Increment &, const oops::Variables &,  Increment &) const {}
    void interpolate(const State &, const oops::Variables &,  State &) const {}
    void smooth(Increment &) const {}
    void save(Increment &) const {}
    void read(const Increment &) {}
    void write(const eckit::Configuration &) const {}

    static const std::string classname() {return "quench::Interpolator";}

 private:
    void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
