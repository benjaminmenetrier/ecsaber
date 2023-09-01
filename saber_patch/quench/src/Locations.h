/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/ObsSpace.h"

namespace quench {

/// Locations class to handle locations for quench model.

class Locations : public util::Printable,
                  private util::ObjectCounter<Locations> {
 public:
  static const std::string classname() {return "quench::Locations";}

  Locations(const ObsSpace & ot,
            const util::DateTime & t1,
            const util::DateTime & t2) {}

  ~Locations() {}

 private:
  void print(std::ostream & os) const {}
};

}  // namespace quench

