/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/range/adaptor/reversed.hpp>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/base/Ensemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelData.h"
#include "oops/interface/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/ECUtilities.h"

namespace saber {

// -----------------------------------------------------------------------------

oops::patch::Variables getActiveVars(const SaberBlockParametersBase &,
                              const oops::patch::Variables &);

// -----------------------------------------------------------------------------

void setMember(eckit::LocalConfiguration & conf,
               const int & member);

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi);

// -----------------------------------------------------------------------------

template<typename MODEL>
eckit::LocalConfiguration readEnsemble(const oops::Geometry<MODEL> & geom,
                                       const oops::patch::Variables & vars,
                                       const oops::State<MODEL> & xb,
                                       const oops::State<MODEL> & fg,
                                       const eckit::LocalConfiguration & inputConf,
                                       const bool & iterativeEnsembleLoading,
                                       std::vector<atlas::FieldSet> & fsetEns) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Prepare ensemble configuration
  oops::Log::info() << "Info     : Prepare ensemble configuration" << std::endl;

  // Create output configuration
  eckit::LocalConfiguration outputConf;

  // Fill output configuration and set ensemble size
  size_t nens = 0;
  size_t ensembleFound = 0;

  // Ensemble of states, perturbation using the mean
  eckit::LocalConfiguration ensembleConf;
  if (inputConf.has("ensemble")) {
    ensembleConf = inputConf.getSubConfiguration("ensemble");
    nens = ensembleConf.getInt("members");
    outputConf.set("ensemble", ensembleConf.getSubConfigurations("state"));
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk
  eckit::LocalConfiguration ensemblePert;
  if (inputConf.has("ensemble pert")) {
    ensemblePert = inputConf.getSubConfiguration("ensemble pert");
    nens = ensemblePert.getInt("members");
    outputConf.set("ensemble", ensemblePert.getSubConfigurations("state"));
    ++ensembleFound;
  }

  // Set ensemble size
  outputConf.set("ensemble size", nens);

  // Check number of ensembles in yaml
  ASSERT(ensembleFound <= 1);

  if (!iterativeEnsembleLoading) {
    // Full ensemble loading
    oops::Log::info() << "Info     : Read full ensemble" << std::endl;

    // Ensemble pointer
    std::unique_ptr<oops::Ensemble<MODEL>> ensemble;

    // Ensemble of states, perturbation using the mean
    if (!ensembleConf.empty()) {
      oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                        << std::endl;
      ensemble.reset(new oops::Ensemble<MODEL>(xb.validTime(), ensembleConf));
      ensemble->linearize(xb, geom);
      for (size_t ie = 0; ie < nens; ++ie) {
        (*ensemble)[ie] *= std::sqrt(static_cast<double>(nens-1));
      }
    }

    // Increment ensemble from increments on disk
    if (!ensemblePert.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;
      ensemble.reset(new oops::Ensemble<MODEL>(xb.validTime(), ensemblePert));
      ensemble->build(xb, geom);
      ensemble->read();
    }

    // Transform Increment into FieldSet
    for (size_t ie = 0; ie < nens; ++ie) {
      atlas::FieldSet fset = util::copyFieldSet((*ensemble)[ie].increment().fieldSet());
      fset.name() = "ensemble member";
      fsetEns.push_back(fset);
    }
  }

  // Return ensemble configuration for iterative ensemble loading
  return outputConf;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readHybridWeight(const oops::Geometry<MODEL> & geom,
                      const oops::patch::Variables & vars,
                      const util::DateTime & date,
                      const eckit::LocalConfiguration & conf,
                      atlas::FieldSet & fset) {
  oops::Log::trace() << "readHybridWeight starting" << std::endl;

  oops::Log::info() << "Info     : Read hybrid weight" << std::endl;

  // Local copy
  eckit::LocalConfiguration localConf(conf);

  // Create Increment
  oops::Increment<MODEL> dx(geom, templatedVars<MODEL>(vars), date);

  // Read file
  dx.read(localConf);

  // Get FieldSet
  fset = util::shareFields(dx.increment().fieldSet());

  oops::Log::trace() << "readHybridWeight done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readEnsembleMember(const oops::Geometry<MODEL> & geom,
                        const oops::patch::Variables & vars,
                        const util::DateTime & date,
                        const eckit::LocalConfiguration & conf,
                        const size_t & ie,
                        atlas::FieldSet & fset) {
  oops::Log::trace() << "readEnsembleMember starting" << std::endl;

  oops::Log::info() << "Info     : Read ensemble member " << ie << std::endl;

  // Fill FieldSet
  size_t ensembleFound = 0;

  if (conf.has("ensemble")) {
    // Ensemble of states passed as increments
    std::vector<eckit::LocalConfiguration> membersConf = conf.getSubConfigurations("ensemble");

    // Read state as increment
    oops::Increment<MODEL> dx(geom, templatedVars<MODEL>(vars), date);
    dx.read(membersConf[ie]);

    // Copy FieldSet
    fset = util::copyFieldSet(dx.increment().fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert")) {
    // Increment ensemble from difference of two states
    std::vector<eckit::LocalConfiguration> membersConf = conf.getSubConfigurations("ensemble");

    // Read Increment
    oops::Increment<MODEL> dx(geom, templatedVars<MODEL>(vars), date);
    dx.read(membersConf[ie]);

    // Get FieldSet
    fset = util::copyFieldSet(dx.increment().fieldSet());

    ++ensembleFound;
  }

  // Check number of ensembles in configuration
  ASSERT(ensembleFound <= 1);

  oops::Log::trace() << "readEnsembleMember done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
