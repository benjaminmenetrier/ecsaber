/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/ConfigFunctions.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {

///  Variables
class Variables : public util::Printable,
                  private util::ObjectCounter<Variables> {
 public:
  static const std::string classname() {return "quench::Variables";}

  explicit Variables(const eckit::Configuration & config) {
    if (util::isVector(config)) {
      vars_ = oops::patch::Variables(config.getStringVector("."));
    } else if (config.has("variables")) {
      vars_ = oops::patch::Variables(config, "variables");
    } else if (config.has("variables list")) {
      vars_ = oops::patch::Variables(config, "variables list");
    } else {
      throw eckit::Exception("wrong variables configuration", Here());
    }
  }
  Variables(const Variables & other) : vars_(other.vars_) {}
  ~Variables() {}

  std::vector<std::string> variablesList() const {return vars_.variables();}

 private:
  void print(std::ostream & os) const {os << vars_ << std::endl;}
  oops::patch::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
