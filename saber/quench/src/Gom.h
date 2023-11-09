/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>

#include "util/abor1_cpp.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class Increment;
  class ObsSpace;
  class State;
  class Variables;

/// Gom class to handle local model values for quench model.

class Gom : public util::Printable,
            private util::ObjectCounter<Gom> {
 public:
  static const std::string classname() {return "quench::Gom";}

  Gom(const ObsSpace &, const Variables &, const Increment &,
        const util::DateTime &, const util::DateTime &)
    {ABORT("not implemented yet");}
  Gom(const ObsSpace &, const Variables &, const State &,
        const util::DateTime &, const util::DateTime &)
    {ABORT("not implemented yet");}

  void zero()
    {ABORT("not implemented yet");}
  void random()
    {ABORT("not implemented yet");}
  double dot_product_with(const Gom &) const {ABORT("not implemented yet"); return 0.0;}
  void read(const eckit::Configuration &)
    {ABORT("not implemented yet");}
  void write(const eckit::Configuration &) const
    {ABORT("not implemented yet");}

 private:
  void print(std::ostream &) const;
};

}  // namespace quench
