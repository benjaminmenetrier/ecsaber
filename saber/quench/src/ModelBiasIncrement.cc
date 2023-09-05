/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/ModelBiasIncrement.h"

#include <iostream>

#include "src/Geometry.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
ModelBiasIncrement::ModelBiasIncrement(const Geometry &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ModelBiasIncrement::ModelBiasIncrement(const ModelBiasIncrement &, const bool) {}
// -----------------------------------------------------------------------------
ModelBiasIncrement::ModelBiasIncrement(const ModelBiasIncrement &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ModelBiasIncrement::~ModelBiasIncrement() {}
// -----------------------------------------------------------------------------

}  // namespace quench
