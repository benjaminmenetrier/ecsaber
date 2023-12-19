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

#include "oops/base/FieldSets.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Variables.h"
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

namespace oops {
  class FieldSet3D;
}

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

void allocateFields(oops::FieldSet3D & fset,
                    const oops::patch::Variables & varsToAllocate,
                    const oops::patch::Variables & varsWithLevels,
                    const atlas::FunctionSpace & functionSpace,
                    const bool haloExchange = true);

// -----------------------------------------------------------------------------

template<typename MODEL>
oops::FieldSets readEnsemble(const oops::Geometry<MODEL> & geom,
                             const oops::patch::Variables & modelvars,
                             const oops::State4D<MODEL> & xb,
                             const oops::State4D<MODEL> & fg,
                             const eckit::LocalConfiguration & inputConf,
                             const bool & iterativeEnsembleLoading,
                             eckit::LocalConfiguration & outputConf) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Prepare ensemble configuration
  oops::Log::info() << "Info     : Prepare ensemble configuration" << std::endl;

  // Fill output configuration and set ensemble size
  size_t nens = 0;
  size_t ensembleFound = 0;

  // Ensemble of states, perturbation using the mean
  std::vector<eckit::LocalConfiguration> ensembleConf;
  if (inputConf.has("ensemble")) {
    if (util::isVector(inputConf.getSubConfiguration("ensemble"))) {
      ensembleConf = inputConf.getSubConfigurations("ensemble");
    } else {
      ensembleConf.push_back(inputConf.getSubConfiguration("ensemble"));
    }
    nens = ensembleConf[0].getInt("members");
    outputConf.set("ensemble", ensembleConf);
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk
  std::vector<eckit::LocalConfiguration> ensemblePert;
  if (inputConf.has("ensemble pert")) {
    if (util::isVector(inputConf.getSubConfiguration("ensemble pert"))) {
      ensemblePert = inputConf.getSubConfigurations("ensemble pert");
    } else {
      ensemblePert.push_back(inputConf.getSubConfiguration("ensemble pert"));
    }
    nens = ensemblePert[0].getInt("members");
    outputConf.set("ensemble", ensemblePert);
    ++ensembleFound;
  }

  // Set ensemble size
  outputConf.set("ensemble size", nens);

  // Check number of ensembles in yaml
  ASSERT(ensembleFound <= 1);

  if (!iterativeEnsembleLoading) {
    // Full ensemble loading
    oops::Log::info() << "Info     : Read full ensemble" << std::endl;

    // Ensemble of states, perturbation using the mean
    if (ensembleConf.size() > 0) {
      oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                        << std::endl;

      for (unsigned jsub = 0; jsub < xb.times().size(); ++jsub) {
        std::shared_ptr<oops::Ensemble<MODEL>> ens_k(new oops::Ensemble<MODEL>(xb[jsub].validTime(),
          ensembleConf[jsub]));
        ens_k->linearize(xb[jsub], geom);
        for (size_t ie = 0; ie < nens; ++ie) {
          (*ens_k)[ie] *= std::sqrt(static_cast<double>(nens-1));
        }
        oops::EnsemblesCollection<MODEL>::getInstance().put(xb[jsub].validTime(), ens_k);
      }
    }

    // Increment ensemble from increments on disk
    if (ensemblePert.size() > 0) {
      oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;

      for (unsigned jsub = 0; jsub < xb.times().size(); ++jsub) {
        std::shared_ptr<oops::Ensemble<MODEL>> ens_k(new oops::Ensemble<MODEL>(xb[jsub].validTime(),
          ensemblePert[jsub]));
        ens_k->build(xb[jsub], geom);
        ens_k->read();
        oops::EnsemblesCollection<MODEL>::getInstance().put(xb[jsub].validTime(), ens_k);
      }
    }

    if (ensembleConf.size() > 0 || ensemblePert.size() > 0) {
      // Transform Ensemble into FieldSets
      oops::Log::info() << "Info     : Transform Ensemble into FieldSets" << std::endl;
      std::vector<int> ensmems(nens);
      std::iota(ensmems.begin(), ensmems.end(), 0);
      oops::FieldSets fsetEns(xb, ensmems);
      return fsetEns;
    }
  }

  // Return empty ensemble if none was returned before
  std::vector<util::DateTime> dates;
  std::vector<int> ensmems;
  oops::FieldSets fsetEns(dates, eckit::mpi::self(), ensmems, eckit::mpi::self());
  return fsetEns;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readHybridWeight(const oops::Geometry<MODEL> & geom,
                      const oops::patch::Variables & vars,
                      const util::DateTime & date,
                      const eckit::LocalConfiguration & conf,
                      oops::FieldSet3D & fset) {
  oops::Log::trace() << "readHybridWeight starting" << std::endl;

  oops::Log::info() << "Info     : Read hybrid weight" << std::endl;

  // Local copy
  eckit::LocalConfiguration localConf(conf);

  // Create variables
  oops::Variables<MODEL> varsT(templatedVarsConf(vars));

  // Create Increment
  oops::Increment<MODEL> dx(geom, varsT, date);

  // Read file
  dx.read(localConf);

  // Get FieldSet
  fset.shallowCopy(dx.increment().fieldSet());

  oops::Log::trace() << "readHybridWeight done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readEnsembleMember(const oops::Geometry<MODEL> & geom,
                        const oops::patch::Variables & vars,
                        const eckit::LocalConfiguration & conf,
                        const size_t & ie,
                        oops::FieldSet3D & fset) {
  oops::Log::trace() << "readEnsembleMember starting" << std::endl;

  oops::Log::info() << "Info     : Read ensemble member " << ie << std::endl;

  // Fill FieldSet
  size_t ensembleFound = 0;

  if (conf.has("ensemble")) {
    // Ensemble of states passed as increments
    std::vector<eckit::LocalConfiguration> ensembleConf = 
      conf.getSubConfigurations("ensemble");
    std::vector<eckit::LocalConfiguration> membersConf =
      ensembleConf[0].getSubConfigurations("state");

    // Create variables
    oops::Variables<MODEL> varsT(templatedVarsConf(vars));

    // Read state as increment
    oops::Increment<MODEL> dx(geom, varsT, fset.validTime());
    dx.read(membersConf[ie]);

    // Copy FieldSet
    fset.deepCopy(dx.increment().fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert")) {
    // Increment ensemble from difference of two states
    std::vector<eckit::LocalConfiguration> ensembleConf
      = conf.getSubConfigurations("ensemble pert");
    std::vector<eckit::LocalConfiguration> membersConf =
      ensembleConf[0].getSubConfigurations("state");

    // Create variables
    oops::Variables<MODEL> varsT(templatedVarsConf(vars));

    // Read Increment
    oops::Increment<MODEL> dx(geom, varsT, fset.validTime());
    dx.read(membersConf[ie]);

    // Get FieldSet
    fset.deepCopy(dx.increment().fieldSet());

    ++ensembleFound;
  }

  // Check number of ensembles in configuration
  ASSERT(ensembleFound <= 1);

  oops::Log::trace() << "readEnsembleMember done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
