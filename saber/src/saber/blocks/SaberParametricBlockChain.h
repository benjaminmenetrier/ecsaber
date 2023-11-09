/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet4D.h"
#include "oops/interface/ModelData.h"
#include "oops/util/ConfigHelpers.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/oops/Utilities.h"

namespace saber {

/// Chain of outer (optional) and not-ensemble central block. Can be used
/// as static error covariance component and as localization for ensemble
/// error covariance.
class SaberParametricBlockChain : public SaberBlockChainBase {
 public:
  template<typename MODEL>
  SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                            const oops::Geometry<MODEL> & dualResGeom,
                            const oops::patch::Variables & outerVars,
                            const oops::FieldSet4D & fset4dXb,
                            const oops::FieldSet4D & fset4dFg,
                            std::vector<atlas::FieldSet> & fsetEns,
                            std::vector<atlas::FieldSet> & fsetDualResEns,
                            const eckit::LocalConfiguration & covarConf,
                            const eckit::Configuration & conf);
  ~SaberParametricBlockChain() = default;

  /// @brief Filter the increment
  void filter(oops::FieldSet4D &) const;

  /// @brief Randomize the increment according to this B matrix.
  void randomize(oops::FieldSet4D &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(oops::FieldSet4D &) const;
  /// @brief Get this B matrix square-root control vector size.
  size_t ctlVecSize() const;
  /// @brief Multiply the control vector by this B matrix square-root.
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const;
  /// @brief Multiply the increment by this B matrix square-root adjoint.
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const;

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::patch::Variables & outerVariables() const {return outerVariables_;}

 private:
  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::patch::Variables outerVariables_;
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  const bool crossTimeCov_;
  std::unique_ptr<SaberCentralBlockBase> centralBlock_;
  const eckit::mpi::Comm & timeComm_;
  size_t size4D_;
  oops::patch::Variables centralVars_;
  atlas::FunctionSpace centralFunctionSpace_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberParametricBlockChain::SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResolutionGeom,
                       const oops::patch::Variables & outerVars,
                       const oops::FieldSet4D & fset4dXb,
                       const oops::FieldSet4D & fset4dFg,
                       // TODO(AS): read inside the block so there is no need to pass
                       // as non-const
                       std::vector<atlas::FieldSet> & fsetEns,
                       std::vector<atlas::FieldSet> & fsetDualResEns,
                       const eckit::LocalConfiguration & covarConf,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.geometry().functionSpace()), outerVariables_(outerVars),
  crossTimeCov_(covarConf.getString("time covariance") == "multivariate duplicated"),
  timeComm_(oops::mpi::myself()), size4D_(fset4dXb.size()) {
  oops::Log::trace() << "SaberParametricBlockChain ctor starting" << std::endl;

  // If needed create outer block chain
  if (conf.has("saber outer blocks")) {
    std::vector<SaberOuterBlockParametersWrapper> cmpOuterBlocksParams;
    for (const auto & cmpOuterBlockConf : conf.getSubConfigurations("saber outer blocks")) {
      SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
      cmpOuterBlockParamsWrapper.deserialize(cmpOuterBlockConf);
      cmpOuterBlocksParams.push_back(cmpOuterBlockParamsWrapper);
    }
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVars,
                          fset4dXb, fset4dFg, fsetEns, covarConf,
                          cmpOuterBlocksParams);
  }

  // Set outer variables and geometry data for central block
  const oops::patch::Variables currentOuterVars = outerBlockChain_ ?
                               outerBlockChain_->innerVars() : outerVariables_;
  const oops::GeometryData geomData(geom.geometry().functionSpace(), geom.geometry().fields(),
    geom.geometry().levelsAreTopDown(), eckit::mpi::comm());
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : geomData;

  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
  saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));

  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

  // Get active variables
  oops::patch::Variables activeVars = getActiveVars(saberCentralBlockParams, currentOuterVars);
  // Check that active variables are present in variables
  for (const auto & var : activeVars.variables()) {
    if (!currentOuterVars.has(var)) {
      throw eckit::UserError("Active variable " + var + " is not present in "
                             "outer variables", Here());
    }
  }

  // Create central block
  centralBlock_ = SaberCentralBlockFactory::create(currentOuterGeom,
                                                   activeVars,
                                                   covarConf,
                                                   saberCentralBlockParams,
                                                   fset4dXb[0],
                                                   fset4dFg[0]);

  // Save central function space and variables
  centralFunctionSpace_ = currentOuterGeom.functionSpace();
  centralVars_ = activeVars;

  // Read and add model fields
  centralBlock_->read(geom, currentOuterVars, fset4dXb[0].validTime());

  // Ensemble configuration
  eckit::LocalConfiguration ensembleConf
         = covarConf.getSubConfiguration("ensemble configuration");
  if (saberCentralBlockParams.doCalibration()) {
    // Block calibration
    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      centralBlock_->iterativeCalibrationInit();

      // Get ensemble size
      size_t nens = ensembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < nens; ++ie) {
        // Read ensemble member
        atlas::FieldSet fset;
        readEnsembleMember(geom,
                           outerVars,
                           fset4dXb[0].validTime(),
                           ensembleConf,
                           ie,
                           fset);

        // Apply outer blocks inverse (all of them)
        oops::Log::info() << "Info     : Apply outer blocks inverse (all of them)" << std::endl;
        if (outerBlockChain_) outerBlockChain_->leftInverseMultiply(fset);

        // Use FieldSet in the central block
        oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
        centralBlock_->iterativeCalibrationUpdate(fset);
      }

      // Finalization
      oops::Log::info() << "Info     : Finalization" << std::endl;
      centralBlock_->iterativeCalibrationFinal();
    } else {
      // Direct calibration
      oops::Log::info() << "Info     : Direct calibration" << std::endl;
      centralBlock_->directCalibration(fsetEns);
    }
  } else if (saberCentralBlockParams.doRead()) {
    // Read data
    oops::Log::info() << "Info     : Read data" << std::endl;
    centralBlock_->read();
  }

  // Dual resolution ensemble
  if (covarConf.has("dual resolution ensemble configuration")) {
    oops::Log::info() << "Info     : Dual resolution setup" << std::endl;

    // Dual resolution setup
    const oops::GeometryData dualResolutionGeomData(geom.geometry().functionSpace(), geom.geometry().fields(),
      geom.geometry().levelsAreTopDown(), eckit::mpi::comm());
    centralBlock_->dualResolutionSetup(dualResolutionGeomData);

    // Ensemble configuration
    eckit::LocalConfiguration dualResolutionEnsembleConf
      = covarConf.getSubConfiguration("dual resolution ensemble configuration");

    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      centralBlock_->iterativeCalibrationInit();

     // Get dual resolution ensemble size
      size_t dualResolutionNens = dualResolutionEnsembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < dualResolutionNens; ++ie) {
        // Read ensemble member
        atlas::FieldSet fset;
        readEnsembleMember(dualResolutionGeom,
                           outerVars,
                           fset4dXb[0].validTime(),
                           dualResolutionEnsembleConf,
                           ie,
                           fset);

        // Use FieldSet in the central block
        oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
        centralBlock_->iterativeCalibrationUpdate(fset);
     }

      // Finalization
      oops::Log::info() << "Info     : Finalization" << std::endl;
      centralBlock_->iterativeCalibrationFinal();
    } else {
      // Direct calibration
      oops::Log::info() << "Info     : Direct calibration" << std::endl;
      centralBlock_->directCalibration(fsetDualResEns);
    }
  }

  // Write calibration data
  if (saberCentralBlockParams.doCalibration()) {
    oops::Log::info() << "Info     : Write calibration data" << std::endl;
    centralBlock_->write(geom, currentOuterVars, fset4dXb[0].validTime());
    centralBlock_->write();
  }

  // Write final ensemble
  if (covarConf.has("output ensemble")) {
    // Get output parameters configuration
    const eckit::LocalConfiguration outputEnsembleConf(covarConf, "output ensemble");

    // Check whether geometry grid is similar to the last outer block inner geometry
    const bool useModelWriter = (util::getGridUid(geom.geometry().functionSpace())
      == util::getGridUid(currentOuterGeom.functionSpace()));

    // Get ensemble size
    size_t ensembleSize = ensembleConf.getInt("ensemble size");

    // Estimate mean
    atlas::FieldSet fsetMean;
    if (iterativeEnsembleLoading) {
      for (size_t ie = 0; ie < ensembleSize; ++ie) {
        // Read member
        atlas::FieldSet fsetMem;
        readEnsembleMember(geom, activeVars, fset4dXb[0].validTime(), ensembleConf, ie, fsetMem);

        // Update mean
        if (ie == 0) {
          fsetMean = util::copyFieldSet(fsetMem);
        } else {
          util::addFieldSets(fsetMean, fsetMem);
        }
      }

      // Normalize mean
      util::multiplyFieldSet(fsetMean, 1.0/static_cast<double>(ensembleSize));
    }

    // Write first member only
    const bool firstMemberOnly = outputEnsembleConf.getBool("first member only", false);
    if (firstMemberOnly) {
      ensembleSize = 1;
    }

    for (size_t ie = 0; ie < ensembleSize; ++ie) {
      oops::Log::info() << "Info     : Write member " << ie << std::endl;

      // Increment pointer
      oops::Increment<MODEL> dx(geom, templatedVars<MODEL>(activeVars), fset4dXb[0].validTime());

      // Get ensemble member
      if (iterativeEnsembleLoading) {
        // Read ensemble member
        readEnsembleMember(geom, activeVars, fset4dXb[0].validTime(), ensembleConf, ie,
                           dx.increment().fieldSet());

        // Remove mean
        util::subtractFieldSets(dx.increment().fieldSet(), fsetMean);

        // Apply outer blocks inverse
        if (outerBlockChain_) outerBlockChain_->leftInverseMultiply(dx.increment().fieldSet());
      } else {
        // Copy member
        dx.increment().fieldSet() = fsetEns[ie];
      }

      // ATLAS fieldset to Increment_
      dx.increment().synchronizeFields();

      if (useModelWriter) {
        // Use model writer

        // Set member index
        eckit::LocalConfiguration outputMemberConf(outputEnsembleConf);
        util::setMember(outputMemberConf, ie+1);

        // Write Increment
        dx.write(outputMemberConf);
        oops::Log::test() << "Norm of ensemble member " << ie << ": " << dx.norm() << std::endl;
      } else {
        // Use generic ATLAS writer
        throw eckit::NotImplemented("generic output ensemble write not implemented yet", Here());
      }
    }
  }

  // Adjoint test
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberCentralBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));

    // Run test
    centralBlock_->adjointTest(currentOuterGeom,
                               activeVars,
                               localAdjointTolerance);
  }

  // Square-root test
  if (covarConf.getBool("square-root test")) {
    // Get tolerance
    const double localSqrtTolerance =
      saberCentralBlockParams.sqrtTolerance.value().get_value_or(
      covarConf.getDouble("square-root tolerance"));

    // Run test
    centralBlock_->sqrtTest(currentOuterGeom,
                            activeVars,
                            localSqrtTolerance);
  }

  oops::Log::trace() << "SaberParametricBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber