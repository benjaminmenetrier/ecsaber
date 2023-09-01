/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/IncrModCtlVec.h"

#include <algorithm>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "util/DateTime.h"
#include "util/Logger.h"

#include "src/ErrorCovariance.h"

using oops::Log;

namespace quench {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrModCtlVec::IncrModCtlVec(const ErrorCovariance &) {}
// -----------------------------------------------------------------------------
IncrModCtlVec::IncrModCtlVec(const ErrorCovariance &,
                             const IncrModCtlVec &) {}
// -----------------------------------------------------------------------------
IncrModCtlVec::IncrModCtlVec(const IncrModCtlVec &,
                             const bool) {}
// -----------------------------------------------------------------------------
IncrModCtlVec::~IncrModCtlVec() {}
// -----------------------------------------------------------------------------
// TODO(Benjamin): for cy46 only, should be removed later
IncrModCtlVec::IncrModCtlVec(const Geometry &,
                             const ErrorCovariance &) {}
IncrModCtlVec::IncrModCtlVec(const Geometry &,
                             const ErrorCovariance &,
                             const IncrModCtlVec &) {}
// -----------------------------------------------------------------------------

}  // namespace quench
