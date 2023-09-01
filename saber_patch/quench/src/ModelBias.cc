/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/ModelBias.h"

#include <iostream>
#include <string>

#include "src/Geometry.h"
#include "src/Model.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const Geometry &, const Model &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const Geometry &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const Geometry &, const ModelBias &) {}
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const ModelBias &, const bool) {}
// -----------------------------------------------------------------------------
ModelBias::~ModelBias() {}
// -----------------------------------------------------------------------------

}  // namespace quench
