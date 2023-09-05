/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#pragma once

#include <omp.h>

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"

#include "util/Logger.h"

namespace saber {

/*! @brief ConvertState application

 * Convert state to another geometry
 */

template <typename MODEL>
class ConvertState : public oops::Application {
  using Geometry_ = oops::Geometry<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using State_ = oops::State<MODEL>;

 public:
  /// Construct a ConvertState application
  ConvertState() {}

  /// Destroy a ConvertState application
  virtual ~ConvertState() {}

  /// Execute a ConvertState application
  int execute(const eckit::Configuration& fullConfig) const {
    // Setup geometry for input and output
    const eckit::LocalConfiguration inputGeomConf(fullConfig, "input geometry");
    const Geometry_ inputGeom(inputGeomConf);
    const eckit::LocalConfiguration outputGeomConf(fullConfig, "output geometry");
    const Geometry_ outputGeom(outputGeomConf);

    // Setup model
    eckit::LocalConfiguration modelConf;
    if (fullConfig.has("model")) {
      modelConf = fullConfig.getSubConfiguration("model");
    }
    const Model_ model(inputGeom, modelConf);

    // List of input and output states
    const std::vector<eckit::LocalConfiguration> statesConf =
      fullConfig.getSubConfigurations("states");
    const int nstates = statesConf.size();

    // Loop over states
    for (int jm = 0; jm < nstates; ++jm) {
      // Print output
      oops::Log::info() << "Converting state " << jm+1 << " of " << nstates << std::endl;

      // Read state
      const eckit::LocalConfiguration inputConf(statesConf[jm], "input");
      State_ xxi(inputGeom, model, inputConf);
      oops::Log::test() << "Input state: " << xxi << std::endl;

      // Copy and change geometry
      State_ xx(outputGeom, xxi);

      // Write state
      const eckit::LocalConfiguration outputConf(statesConf[jm], "output");
      xx.write(outputConf);

      oops::Log::test() << "Output state: " << xx << std::endl;
    }

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::ConvertState<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace saber
