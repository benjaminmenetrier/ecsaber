/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

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
      config_.set("variables", eckit::LocalConfiguration(config));
    } else {
      config_ = eckit::LocalConfiguration(config);
    }
    vars_ = oops::patch::Variables(config_, "variables");
  }
  Variables(const Variables & other) : config_(other.config_), vars_(other.vars_) {}
  ~Variables() {};

  const eckit::LocalConfiguration & config() const {return config_;}
  std::vector<std::string> varlist() const {return vars_.variables();}

 private:
  void print(std::ostream & os) const {os << vars_ << std::endl;}
  eckit::LocalConfiguration config_;
  oops::patch::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
