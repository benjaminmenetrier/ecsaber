/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <omp.h>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "oops/runs/Application.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"
#include "saber/oops/ECUtilities.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ProcessPerts application.
class ProcessPertsParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ProcessPertsParameters, oops::Parameters)

 public:
  /// Geometry parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"resolution", this};

  /// Background parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> background{"Background", this};

  oops::RequiredParameter<util::DateTime> date{"date", this};
  oops::RequiredParameter<eckit::LocalConfiguration> inputVariables{"input variables", this};

  /// Background error covariance model.
  oops::RequiredParameter<eckit::LocalConfiguration> filterCovarianceBlockConf{"saber filter blocks", this};

  /// Where to read input ensemble: From states or perturbations
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};

  /// Where to write low-pass filtered perturbations
  oops::OptionalParameter<eckit::LocalConfiguration>
    lowpassPerturbations{"low pass perturbations", this};

  /// Where to write high-pass filtered perturbations
  oops::RequiredParameter<eckit::LocalConfiguration>
    outputPerturbations{"output perturbations", this};
};

// -----------------------------------------------------------------------------

/*! @brief Error covariance training application

 * Use input fields and ensemble data to perform the calibration of
 * a background error covariance matrix
 */

template <typename MODEL>
class ProcessPerts : public oops::Application {
  using CovarianceFactory_ = oops::CovarianceFactory<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using Localization_ = oops::Localization<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using CovarianceBase_ = oops::ModelSpaceCovarianceBase<MODEL>;
  using State_ = oops::State<MODEL>;
  using State4D_ = oops::State4D<MODEL>;
  using Variables_ = oops::Variables<MODEL>;

 public:
// -----------------------------------------------------------------------------
  ProcessPerts() {
    instantiateCovarFactory<MODEL>();
    instantiateLocalizationFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~ProcessPerts() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    // Deserialize parameters
    ProcessPertsParameters params;
    params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Define space and time communicators
    const eckit::mpi::Comm * commSpace = &eckit::mpi::comm();

    // Get number of MPI tasks and OpenMP threads
    size_t ntasks = commSpace->size();
    size_t nthreads = 1;
#ifdef _OPENMP
    # pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }
#endif

    // Replace patterns in full configuration and deserialize parameters
    eckit::LocalConfiguration fullConfigUpdated(fullConfig);
    util::seekAndReplace(fullConfigUpdated, "_MPI_", std::to_string(ntasks));
    util::seekAndReplace(fullConfigUpdated, "_OMP_", std::to_string(nthreads));
    params.deserialize(fullConfigUpdated);

    // Set precision for test channel
    oops::Log::test() << std::scientific
                      << std::setprecision(std::numeric_limits<double>::digits10+1);

    // Setup geometry
    const Geometry_ geom(params.geometry);

    // Setup model
    eckit::LocalConfiguration modelConf;
    if (fullConfigUpdated.has("model")) {
      modelConf = fullConfigUpdated.getSubConfiguration("model");
    }
    const Model_ model(geom, modelConf);

    // Setup background state
    const State4D_ xx(params.background, geom, model);
    oops::FieldSet4D fsetXb(xx);
    oops::FieldSet4D fsetFg(xx);

    // Setup variables
    const Variables_ statevars(params.inputVariables.value());

    // Setup time
    const util::DateTime time = xx[0].validTime();

    eckit::LocalConfiguration filterCovarianceBlockConf(fullConfig, "saber filter blocks");
    std::unique_ptr<SaberParametricBlockChain> saberFilterBlocks;

    // List of output increments
    const auto & lowpassPerturbations = params.lowpassPerturbations.value();
    const auto & incrementsWriteParams =
      params.outputPerturbations.value();

    Variables_ incVarsT(params.inputVariables.value());
    oops::patch::Variables incVars(incVarsT.variables().variablesList());
    // Initialize outer variables
    const std::vector<std::size_t> vlevs = geom.geometry().variableSizes(incVarsT.variables());
    for (std::size_t i = 0; i < vlevs.size() ; ++i) {
      incVars.addMetaData(incVars[i], "levels", vlevs[i]);
    }

    std::vector<atlas::FieldSet> fsetEns;
    std::vector<atlas::FieldSet> dualResFsetEns;
    eckit::LocalConfiguration covarConf;
    covarConf.set("iterative ensemble loading", false);
    covarConf.set("inverse test", false);
    covarConf.set("adjoint test", false);
    covarConf.set("square-root test", false);
    covarConf.set("covariance model", "SABER");
    covarConf.set("time covariance", "");

    // Initialize filter blockchain
    saberFilterBlocks = std::make_unique<SaberParametricBlockChain>(geom, geom,
      incVars, fsetXb, fsetFg, fsetEns, dualResFsetEns,
      covarConf, filterCovarianceBlockConf);

    // Yaml validation
    // TODO(Mayeul): Move this do an override of deserialize
    if (((params.ensemble.value() == boost::none) &&
        (params.ensemblePert.value() == boost::none)) ||
        ((params.ensemble.value() != boost::none) &&
        (params.ensemblePert.value() != boost::none)))
    {
      throw eckit::UserError(
       "Require either input states or input perturbations to be set in yaml",
       Here());
    }

    // Read input ensemble
    const bool iterativeEnsembleLoading = false;
    std::vector<atlas::FieldSet> fsetEnsI;
    eckit::LocalConfiguration ensembleConf(fullConfig);
    readEnsemble<MODEL>(geom,
                        incVars,
                        xx[0], xx[0],
                        ensembleConf,
                        iterativeEnsembleLoading,
                        fsetEnsI);
    int nincrements = fsetEnsI.size();

    //  Loop over perturbations
    for (int jm = 0; jm < nincrements; ++jm) {
      //  Read ensemble member perturbation
      atlas::FieldSet fsetI = fsetEnsI[jm];
      Increment_ dxI(geom, incVarsT, time);
      dxI.zero();
      dxI.increment().fieldSet() = fsetI;
      dxI.increment().synchronizeFields();

      //  Copy perturbation
      atlas::FieldSet fset;
      util::copyFieldSet(fsetI, fset);

      oops::Log::test() << "Norm of perturbation : member  " << jm+1
                        << ": " << dxI.norm() << std::endl;

      oops::FieldSet4D fset4dDxI(oops::FieldSet3D{fsetI, time, eckit::mpi::comm()});
      oops::FieldSet4D fset4dDx(oops::FieldSet3D{fset, time, eckit::mpi::comm()});

      // Apply filter blocks
      saberFilterBlocks->filter(fset4dDx);

      if (lowpassPerturbations != boost::none) {
        Increment_ dxLowPass(geom, incVarsT, time);
        dxLowPass.zero();
        dxLowPass.increment().fieldSet() = fset4dDx[0].fieldSet();
        dxLowPass.increment().synchronizeFields();

        auto lowpassPerturbationsUpdated(*lowpassPerturbations);
        setMember(lowpassPerturbationsUpdated,jm+1);
        dxLowPass.write(lowpassPerturbationsUpdated);
        oops::Log::test() << "Norm of low pass perturbation : member  " << jm+1
                          << ": " << dxLowPass.norm() << std::endl;
      }

      // High pass = full pert - low pass
      fset4dDx *= -1.0;
      fset4dDxI += fset4dDx;

      // Write high pass
      dxI.increment().fieldSet() = fset4dDxI[0].fieldSet();
      dxI.increment().synchronizeFields();

      auto incrementsWriteParamsUpdated(incrementsWriteParams);
      setMember(incrementsWriteParamsUpdated,jm+1);
      dxI.write(incrementsWriteParamsUpdated);

      oops::Log::test() << "Norm of high pass perturbation : member  " << jm+1
                        << ": " << dxI.norm() << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ProcessPerts<" + MODEL::name() + ">";
  }
};

}  // namespace saber
