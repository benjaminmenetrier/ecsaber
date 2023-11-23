/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/ECUtilities.h"

namespace saber {

// -----------------------------------------------------------------------------

eckit::LocalConfiguration templatedVarsConf(const oops::patch::Variables & vars) {
  eckit::LocalConfiguration varConf;
  varConf.set("variables list", vars.variables());
  return varConf;
}

// -----------------------------------------------------------------------------

}  // namespace saber
