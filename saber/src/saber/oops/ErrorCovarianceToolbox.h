/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#pragma once

#include <omp.h>

#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"

namespace saber {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ErrorCovarianceToolbox application.
template <typename MODEL> class ErrorCovarianceToolboxParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceToolboxParameters, oops::Parameters)

 public:
  /// Geometry parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"resolution", this};

  /// Background parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> background{"Background", this};

  /// Background error covariance model.
  oops::RequiredParameter<eckit::LocalConfiguration> backgroundError{"Covariance", this};

  /// Geometry parameters.
  oops::Parameter<bool> parallel{"parallel subwindows", true, this};

  /// Outer variables parameters
  oops::OptionalParameter<oops::Variables> incrementVars{"increment variables", this};

  /// Dirac location/variables parameters.
  oops::OptionalParameter<eckit::LocalConfiguration> dirac{"dirac", this};

  /// Diagnostic location/variables parameters.
  oops::OptionalParameter<eckit::LocalConfiguration> diagnostic{"diagnostic points", this};

  /// Where to write the output(s) of Dirac tests
  oops::OptionalParameter<eckit::LocalConfiguration> outputDirac{"output dirac", this};

  /// Whether and where to write the perturbations generated by the randomization.
  oops::OptionalParameter<eckit::LocalConfiguration> outputPerturbations{"output perturbations",
    this};

  /// Whether and where to write the states generated by the randomization.
  oops::OptionalParameter<eckit::LocalConfiguration> outputStates{"output states", this};

  /// Where to write the output of randomized variance.
  oops::OptionalParameter<eckit::LocalConfiguration> outputVariance{"output variance", this};
};

// -----------------------------------------------------------------------------

/*! @brief Error covariance training application

 * Use input fields and ensemble data to perform the calibration of
 * a background error covariance matrix
 */

template <typename MODEL>
class ErrorCovarianceToolbox : public oops::Application {
  using CovarianceFactory_ = oops::CovarianceFactory<MODEL>;
  using Geometry_ = oops::Geometry<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using Localization_ = oops::Localization<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using CovarianceBase_ = oops::ModelSpaceCovarianceBase<MODEL>;
  using State_ = oops::State<MODEL>;
  using State4D_ = oops::State4D<MODEL>;
  using ErrorCovarianceToolboxParameters_ = ErrorCovarianceToolboxParameters<MODEL>;

 public:
// -----------------------------------------------------------------------------
  ErrorCovarianceToolbox() {
    instantiateCovarFactory<MODEL>();
    instantiateLocalizationFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~ErrorCovarianceToolbox() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    // Deserialize parameters
    ErrorCovarianceToolboxParameters_ params;
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

    // Setup variables
    ASSERT(params.incrementVars.value() != boost::none);
    oops::Variables vars(*params.incrementVars.value());

    // Setup time
    util::DateTime time = xx[0].validTime();

    // Background error covariance parameters
    const eckit::LocalConfiguration & covarParams
      = params.backgroundError.value();

    // Dirac test
    const auto & diracParams = params.dirac.value();
    if (diracParams != boost::none) {
      // Setup Dirac field
      Increment4D_ dxi(geom, vars, xx.times());
      dxi.dirac(*diracParams);
      oops::Log::test() << "Input Dirac increment:" << dxi << std::endl;

      // Test configuration
      eckit::LocalConfiguration testConf;

      // Add dirac configuration
      testConf.set("dirac", *diracParams);

      // Add diagnostic print configuration
      const auto & diagnostic = params.diagnostic.value();
      if (diagnostic != boost::none) {
        testConf.set("diagnostic points", *diagnostic);
      }

      // Dirac output parameters
      const auto & outputDirac = params.outputDirac.value();

      // Update parameters
      auto outputDiracUpdated(*outputDirac);
      setMPI(outputDiracUpdated, ntasks);

      // Add output Dirac configuration
      eckit::LocalConfiguration outputConf(outputDiracUpdated);
      testConf.set("output dirac", outputConf);

      // Apply B matrix components recursively
      std::string id;
      dirac(covarParams, testConf, id, geom, vars, xx, dxi);
    }

    const auto & randomizationSize = covarParams.getInt("randomization size", 0);
    if ((diracParams == boost::none) || (randomizationSize > 0)) {
      // Background error covariance training
      std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                            covarParams, geom, vars, xx[0]));

      // Linearize
      eckit::LocalConfiguration linConf;
      const std::string covarianceModel(covarParams.getString("covariance"));
      if (covarianceModel == "hybrid") {
        eckit::LocalConfiguration jbConf;
        jbConf.set("Covariance", covarParams);
        linConf.set("Jb", jbConf);
      } else if (covarianceModel == "ensemble") {
        linConf.set("ensemble_covariance", covarParams);
      } else {
        linConf = covarParams;
      }
      Bmat->linearize(xx[0], geom, linConf);

      // Randomization
      randomization(params, geom, vars, xx, Bmat, ntasks);
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ErrorCovarianceToolbox<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
// The passed geometry should be consistent with the passed increment
  void print_value_at_positions(const eckit::LocalConfiguration & diagConf,
                                const Geometry_ & geom,
                                const Increment4D_ & data) const {
    oops::Log::trace() << appname() << "::print_value_at_position starting" << std::endl;

    // Create diagnostic field
    Increment4D_ diagPoints(data);
    diagPoints.dirac(diagConf);

    // Get diagnostic values
    for (size_t jj = 0; jj < data.size(); ++jj) {
      util::printDiagValues(diagPoints.commTime(),
                            geom.getComm(),
                            geom.functionSpace(),
                            data[jj].fieldSet(),
                            diagPoints[jj].fieldSet());
    }

    oops::Log::trace() << appname() << "::print_value_at_position done" << std::endl;
  }
// -----------------------------------------------------------------------------
// The passed geometry/variables should be consistent with the passed increment
  void dirac(const eckit::LocalConfiguration & covarConf,
             const eckit::LocalConfiguration & testConf,
             std::string & id,
             const Geometry_ & geom,
             const oops::Variables & vars,
             const State4D_ & xx,
             const Increment4D_ & dxi) const {

    // Define output increment
    Increment4D_ dxo(dxi, false);

    // Covariance
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                          covarConf, geom, vars, xx[0]));

    // Linearize
    eckit::LocalConfiguration linConf;
    const std::string covarianceModel(covarConf.getString("covariance"));
    if (covarianceModel == "hybrid") {
      eckit::LocalConfiguration jbConf;
      jbConf.set("Covariance", covarConf);
      linConf.set("Jb", jbConf);
    } else if (covarianceModel == "ensemble") {
      linConf.set("ensemble_covariance", covarConf);
    } else {
      linConf = covarConf;
    }
    Bmat->linearize(xx[0], geom, linConf);

    // Multiply
    Bmat->multiply(dxi[0], dxo[0]);

    // Update ID
    if (id != "") id.append("_");
    id.append(Bmat->covarianceModel());

    if (testConf.has("diagnostic points")) {
      oops::Log::test() << "Covariance(" << id << ") diagnostics:" << std::endl;

      // Print variances
      oops::Log::test() << "- Variances at Dirac points:" << std::endl;
      print_value_at_positions(testConf.getSubConfiguration("dirac"), geom, dxo);

      // Print covariances
      oops::Log::test() << "- Covariances at diagnostic points:" << std::endl;
      print_value_at_positions(testConf.getSubConfiguration("diagnostic points"), geom, dxo);
    }

    // Copy configuration
    eckit::LocalConfiguration outputBConf(testConf.getSubConfiguration("output dirac"));

    // Seek and replace %id% with id, recursively
    util::seekAndReplace(outputBConf, "%id%", id);

    // Write output increment
    dxo[0].write(outputBConf);
    oops::Log::test() << "Covariance(" << id << ") * Increment:" << dxo << std::endl;

std::cout << "toto 1 " << id << std::endl;
    // Look for hybrid or ensemble covariance models
    if (covarianceModel == "hybrid") {
      eckit::LocalConfiguration staticConfig(covarConf, "static_covariance");
      std::string staticID = "hybrid1";
      dirac(staticConfig, testConf, staticID, geom, vars, xx, dxi);
      eckit::LocalConfiguration ensembleConfig(covarConf, "ensemble_covariance");
      std::string ensembleID = "hybrid2";
      dirac(ensembleConfig, testConf, ensembleID, geom, vars, xx, dxi);
    }
    if (covarianceModel == "ensemble" && covarConf.has("localization")) {
      // Localization configuration
      eckit::LocalConfiguration locConfig(covarConf.getSubConfiguration("localization"));
      locConfig.set("date", xx[0].validTime().toString());
      locConfig.set("variables", vars.variables());

      // Define output increment
      Increment4D_ dxo(dxi);

      // Setup localization
      Localization_ Lmat(geom, locConfig);

      // Apply localization
      Lmat.multiply(dxo[0]);

      // Update ID
      std::string idL(id);
      idL.append("_localization");

      // Print localization
      oops::Log::test() << "Localization(" << idL << ") diagnostics:" << std::endl;
      oops::Log::test() << "- Localization at zero separation:" << std::endl;
      print_value_at_positions(testConf.getSubConfiguration("dirac"), geom, dxo);
      if (testConf.has("diagnostic points")) {
        oops::Log::test() << "- Localization at diagnostic points:" << std::endl;
        print_value_at_positions(testConf.getSubConfiguration("diagnostic points"), geom, dxo);
      }

      // Copy configuration
      eckit::LocalConfiguration outputLConf(testConf.getSubConfiguration("output dirac"));

      // Seek and replace %id% with id, recursively
      util::seekAndReplace(outputLConf, "%id%", idL);

      // Write output increment
      dxo[0].write(outputLConf);
      oops::Log::test() << "Localization(" << id << ") * Increment:" << dxo << std::endl;
    }
std::cout << "toto 2 " << id << std::endl;
  }
// -----------------------------------------------------------------------------
  void randomization(const ErrorCovarianceToolboxParameters_ & params,
                     const Geometry_ & geom,
                     const oops::Variables & vars,
                     const State4D_ & xx,
                     const std::unique_ptr<CovarianceBase_> & Bmat,
                     const size_t & ntasks) const {
    if (Bmat->randomizationSize() > 0) {
      oops::Log::info() << "Info     : " << std::endl;
      oops::Log::info() << "Info     : Generate perturbations:" << std::endl;
      oops::Log::info() << "Info     : -----------------------" << std::endl;

      // Create increments
      Increment4D_ dx(geom, vars, xx.times());
      Increment4D_ dxsq(geom, vars, xx.times());
      Increment4D_ variance(geom, vars, xx.times());

      // Initialize variance
      variance.zero();

      // Create empty ensemble
      std::vector<Increment_> ens;

      // Output options
      const auto & outputPerturbations = params.outputPerturbations.value();
      const auto & outputStates = params.outputStates.value();
      const auto & outputVariance = params.outputVariance.value();

      for (size_t jm = 0; jm < Bmat->randomizationSize(); ++jm) {
        // Generate member
        oops::Log::info() << "Info     : Member " << jm << std::endl;
        Bmat->randomize(dx[0]);

        if ((outputPerturbations != boost::none) || (outputStates != boost::none)) {
          // Save member
          ens.push_back(dx[0]);
        }

        // Square perturbation
        dxsq = dx;
        dxsq.schur_product_with(dx);

        // Update variance
        variance += dxsq;
      }
      oops::Log::info() << "Info     : " << std::endl;

      if ((outputPerturbations != boost::none) || (outputStates != boost::none)) {
        oops::Log::info() << "Info     : Write states and/or perturbations:" << std::endl;
        oops::Log::info() << "Info     : ----------------------------------" << std::endl;
        oops::Log::info() << "Info     : " << std::endl;
        for (size_t jm = 0; jm < Bmat->randomizationSize(); ++jm) {
          oops::Log::test() << "Member " << jm << ": " << ens[jm] << std::endl;

          if (outputPerturbations != boost::none) {
            // Update parameters
            auto outputPerturbationsUpdated(*outputPerturbations);
            setMember(outputPerturbationsUpdated, jm+1);
            setMPI(outputPerturbationsUpdated, ntasks);

            // Write perturbation
            ens[jm].write(outputPerturbationsUpdated);
          }

          if (outputStates != boost::none) {
            // Update parameters
            auto outputStatesUpdated(*outputStates);
            setMember(outputStatesUpdated, jm+1);
            setMPI(outputStatesUpdated, ntasks);

            // Add background state to perturbation
            State_ xp(xx[0]);
            xp += ens[jm];

            // Write state
            xp.write(outputStatesUpdated);
          }

          oops::Log::info() << "Info     : " << std::endl;
        }
      }

      if (outputVariance != boost::none) {
        oops::Log::info() << "Info     : Write randomized variance:" << std::endl;
        oops::Log::info() << "Info     : --------------------------" << std::endl;
        oops::Log::info() << "Info     : " << std::endl;
        if (Bmat->randomizationSize() > 1) {
          // Normalize variance
          double rk_norm = 1.0/static_cast<double>(Bmat->randomizationSize());
          variance *= rk_norm;
        }

        // Update parameters
        auto outputVarianceUpdated = *outputVariance;
        setMPI(outputVarianceUpdated, ntasks);

        // Write variance
        variance[0].write(outputVarianceUpdated);
        oops::Log::test() << "Randomized variance: " << variance << std::endl;
      }
    }
  }
// -----------------------------------------------------------------------------
};

}  // namespace saber
