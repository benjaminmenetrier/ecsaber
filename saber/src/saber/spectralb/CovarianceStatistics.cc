/* 
 * (C) Crown Copyright 2017-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0 
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "saber/spectralb/CovarianceStatistics.h"

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "saber/spectralb/CovarianceStatisticsUtils.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

CovStat_ErrorCov::CovStat_ErrorCov(const oops::patch::Variables & vars,
                                   const Parameters_ & params) :
  covarianceFileName_(params.covarianceFile),
  modelLevels_(vars.getLevels(vars[0])),
  nSpectralBinsFull_(getNSpectralBinsFull(params)),
  spectralUMatrices_(createUMatrices(vars, modelLevels_,
                                     nSpectralBinsFull_, params)),
  spectralVerticalCovariances_(createSpectralCovariances(
                               vars, modelLevels_, nSpectralBinsFull_,
                               spectralUMatrices_, params)),
  verticalSD_(createVerticalSD(vars,
                               spectralVerticalCovariances_)),
  spectralCorrelUMatrices_(createCorrelUMatrices(vars,
                                                 spectralVerticalCovariances_,
                                                 spectralUMatrices_,
                                                 verticalSD_)),
  spectralVerticalCorrelations_(createSpectralCorrelations(
                                vars, spectralVerticalCovariances_,
                                verticalSD_))
{
  for (const std::string & s : vars.variables()) {
    if (static_cast<int>(vars.getLevels(s)) != modelLevels_) {
      throw eckit::UnexpectedState("spectral covariance block assumes all fields have "
                                   "same number of model levels");
    }
  }
}

void CovStat_ErrorCov::print(std::ostream & os) const {
  oops::Log::trace() <<
    "Covariance Statistics (SABER spectral B, error covariance) print starting" << std::endl;

  os << std::endl << "  covstats print";
  oops::Log::trace() <<
    "Covariance Statistics (SABER spectral B, error covariance) print done" << std::endl;
}

}  // namespace spectralb
}  // namespace saber
