/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/base/Variables.h"

namespace saber {

// -----------------------------------------------------------------------------

eckit::LocalConfiguration templatedVarsConf(const oops::patch::Variables &);

// -----------------------------------------------------------------------------

template<typename MODEL>
void dirac4D(const eckit::Configuration & conf,
             oops::Increment4D<MODEL> & incr4D) {
  if (incr4D.first() == incr4D.last()) {
    incr4D[0].increment().dirac(conf);
  } else {
    const std::vector<eckit::LocalConfiguration> confs = conf.getSubConfigurations();
    ASSERT(incr4D.last()-incr4D.first()+1 == confs.size());
    for (int jt = incr4D.first(); jt <= incr4D.last(); ++jt) {
      if (!confs[jt].empty()) {
        incr4D[jt].increment().dirac(confs[jt]);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
