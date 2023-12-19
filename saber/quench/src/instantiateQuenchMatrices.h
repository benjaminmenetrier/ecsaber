/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include "oops/interface/GenericMatrix.h"

#include "src/HybridWeight.h"
#include "src/Traits.h"

namespace quench {
// -----------------------------------------------------------------------------

void instantiateQuenchMatrices() {
  static oops::GenericMatrixMaker<quench::Traits, HybridWeight> makerEnsWgt_("hybrid_weight");
}

// -----------------------------------------------------------------------------
}  // namespace quench
