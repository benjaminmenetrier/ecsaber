/*
 * (C) Copyright 2021 UCAR
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/assimilation/IncrCtlVec.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberEnsembleBlockChain.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

#include "util/abor1_cpp.h"
#include "util/Logger.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL> {
  using Geometry_ = oops::Geometry<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Increment4D_ = oops::Increment4D<MODEL>;
  using IncrCtlVec_ = oops::IncrCtlVec<MODEL>;
  using IncrEnsCtlVec_ = oops::IncrEnsCtlVec<MODEL>;
  using IncrModCtlVec_ = oops::IncrModCtlVec<MODEL>;
  using State_ = oops::State<MODEL>;
  using State4D_ = oops::State4D<MODEL>;

 public:
  typedef ErrorCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &, const eckit::Configuration &,
                  const State_ &);
  ~ErrorCovariance();

  // Methods
  void linearize(const State_ &, const Geometry_ &) override;
  void multiply(const Increment_ &, Increment_ &) const override;
  void inverseMultiply(const Increment_ &, Increment_ &) const override;
  void multiplySqrt(const IncrCtlVec_ &, Increment_ &) const override;
  void multiplySqrtTrans(const Increment_ &, IncrCtlVec_ &) const override;
  void randomize(Increment_ &) const override;

  // Control Vector
  IncrModCtlVec_ *newIncrModCtlVec() const override {
    return new IncrModCtlVec_(this->ctlVecSize());
  }
  IncrEnsCtlVec_ *newIncrEnsCtlVec() const override {
    return new IncrEnsCtlVec_();
  }

 private:
  ErrorCovariance(const ErrorCovariance&);
  ErrorCovariance& operator=(const ErrorCovariance&);

  size_t ctlVecSize() const;

  void print(std::ostream &) const;

  /// Chain of outer blocks applied to all components of hybrid covariances.
  /// Not initialized for non-hybrid covariances.
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// Vector of hybrid B components (one element for non-hybrid case).
  std::vector<std::unique_ptr<SaberBlockChainBase>> hybridBlockChain_;
  /// Vector of scalar weights for hybrid B components (one element, equal to
  /// 1.0 for non-hybrid case).
  std::vector<double> hybridScalarWeightSqrt_;
  /// Vector of field weights for hybrid B components (one element, empty
  /// fieldset for non-hybrid case).
  std::vector<atlas::FieldSet> hybridFieldWeightSqrt_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const eckit::Configuration & config,
                                        const State_ & xb3D)
  : oops::ModelSpaceCovarianceBase<MODEL>(config)
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // JEDI compatibility
  State4D_ xb;
  xb.push_back(xb3D);
  State4D_ fg(xb);
  Parameters_ params;
  params.validateAndDeserialize(config);

  // Local copy of background and first guess that can undergo interpolation
  oops::FieldSet4D fset4dXb(xb);
  oops::FieldSet4D fset4dFg(fg);

  // Extend background and first guess with geometry fields
  // TODO(Benjamin, Marek, Mayeul, ?)

  // Initialize outer variables
  const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
  oops::Variables outerVars(incVars.variables());
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    outerVars.addMetaData(outerVars[i], "levels", vlevs[i]);
  }

  // Create covariance configuration
  eckit::LocalConfiguration covarConf;
  covarConf.set("adjoint test", params.adjointTest.value());
  covarConf.set("adjoint tolerance", params.adjointTolerance.value());
  covarConf.set("inverse test", params.inverseTest.value());
  covarConf.set("inverse tolerance", params.inverseTolerance.value());
  covarConf.set("square-root test", params.sqrtTest.value());
  covarConf.set("square-root tolerance", params.sqrtTolerance.value());
  covarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());
  covarConf.set("time covariance", params.timeCovariance.value());

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = params.iterativeEnsembleLoading.value();

  // Initialize ensembles as vector of FieldSets
  std::vector<atlas::FieldSet> fsetEns;
  // Read ensemble (for non-iterative ensemble loading)
  eckit::LocalConfiguration ensembleConf = readEnsemble(geom,
                                                        outerVars,
                                                        xb[0],
                                                        fg[0],
                                                        params.toConfiguration(),
                                                        iterativeEnsembleLoading,
                                                        fsetEns);
  covarConf.set("ensemble configuration", ensembleConf);
  // Read dual resolution ensemble if needed
  const auto & dualResParams = params.dualResParams.value();
  const oops::Geometry<MODEL> * dualResGeom = &geom;
  std::vector<atlas::FieldSet> fsetDualResEns;
  if (dualResParams != boost::none) {
    const auto & dualResGeomConf = dualResParams->geometry.value();
    if (dualResGeomConf != boost::none) {
      // Create dualRes geometry
      dualResGeom = new oops::Geometry<MODEL>(*dualResGeomConf);
    }
    // Background and first guess at dual resolution geometry
    oops::State<MODEL> xbDualRes(*dualResGeom, xb[0]);
    oops::State<MODEL> fgDualRes(*dualResGeom, fg[0]);
    // Read dual resolution ensemble
    eckit::LocalConfiguration dualResEnsembleConf
      = readEnsemble(*dualResGeom,
                     outerVars,
                     xbDualRes,
                     fgDualRes,
                     dualResParams->toConfiguration(),
                     iterativeEnsembleLoading,
                     fsetDualResEns);

    // Add dual resolution ensemble configuration
    covarConf.set("dual resolution ensemble configuration", dualResEnsembleConf);
  }

  // Add ensemble output
  const auto & outputEnsemble = params.outputEnsemble.value();
  if (outputEnsemble != boost::none) {
    covarConf.set("output ensemble", *outputEnsemble);
  }

  const SaberBlockParametersBase & saberCentralBlockParams =
    params.saberCentralBlockParams.value().saberCentralBlockParameters;
  // Build covariance blocks: hybrid covariance case
  if (saberCentralBlockParams.saberBlockName.value() == "Hybrid") {
    // Build common (for all hybrid components) outer blocks if they exist
    const auto & saberOuterBlocksParams = params.saberOuterBlocksParams.value();
    if (saberOuterBlocksParams != boost::none) {
      outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom,
                       outerVars,
                       fset4dXb,
                       fset4dFg,
                       fsetEns,
                       covarConf,
                       *saberOuterBlocksParams);
      outerVars = outerBlockChain_->innerVars();
    }

    // Hybrid central block
    eckit::LocalConfiguration hybridConf = saberCentralBlockParams.toConfiguration();

    // Create block geometry (needed for ensemble reading)
    const oops::Geometry<MODEL> * hybridGeom = &geom;
    if (hybridConf.has("geometry")) {
      hybridGeom = new oops::Geometry<MODEL>(hybridConf.getSubConfiguration("geometry"));
    }

    // Loop over components
    for (const auto & cmp : hybridConf.getSubConfigurations("components")) {
      // Initialize component outer variables
      // TODO(AS): this should be either outerVars or outerBlockChain_->innerVars();
      oops::Variables cmpOuterVars(outerVars);

      // Set weight
      eckit::LocalConfiguration weightConf = cmp.getSubConfiguration("weight");
      // Scalar weight
      hybridScalarWeightSqrt_.push_back(std::sqrt(weightConf.getDouble("value", 1.0)));
      // File-base weight
      atlas::FieldSet fsetWeight;
      if (weightConf.has("file")) {
        // File-base weight
        readHybridWeight(*hybridGeom,
                         outerVars,
                         xb[0].validTime(),
                         weightConf.getSubConfiguration("file"),
                         fsetWeight);
        util::sqrtFieldSet(fsetWeight);
      }
      hybridFieldWeightSqrt_.push_back(fsetWeight);

      // Set covariance
      eckit::LocalConfiguration cmpConf = cmp.getSubConfiguration("covariance");

      // Initialize ensembles as vector of FieldSets
      std::vector<atlas::FieldSet> fset4dCmpEns;
      // Read ensemble
      eckit::LocalConfiguration cmpEnsembleConf
         = readEnsemble(*hybridGeom,
                        cmpOuterVars,
                        xb[0],
                        fg[0],
                        cmpConf,
                        params.iterativeEnsembleLoading.value(),
                        fset4dCmpEns);

      // Create internal configuration
      eckit::LocalConfiguration cmpCovarConf;
      cmpCovarConf.set("ensemble configuration", cmpEnsembleConf);
      cmpCovarConf.set("adjoint test", params.adjointTest.value());
      cmpCovarConf.set("adjoint tolerance", params.adjointTolerance.value());
      cmpCovarConf.set("inverse test", params.inverseTest.value());
      cmpCovarConf.set("inverse tolerance", params.inverseTolerance.value());
      cmpCovarConf.set("square-root test", params.sqrtTest.value());
      cmpCovarConf.set("square-root tolerance", params.sqrtTolerance.value());
      cmpCovarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());
      cmpCovarConf.set("time covariance", params.timeCovariance.value());

      SaberCentralBlockParametersWrapper cmpCentralBlockParamsWrapper;
      cmpCentralBlockParamsWrapper.deserialize(cmpConf.getSubConfiguration("saber central block"));
      const auto & centralBlockParams =
                   cmpCentralBlockParamsWrapper.saberCentralBlockParameters.value();
      // TODO(AS): move construction of BlockChain to factory method or function
      if (centralBlockParams.saberBlockName.value() == "Ensemble") {
        hybridBlockChain_.push_back(std::make_unique<SaberEnsembleBlockChain>(*hybridGeom,
                          *dualResGeom,
                          cmpOuterVars,
                          fset4dXb,
                          fset4dFg,
                          fset4dCmpEns,
                          fsetDualResEns,
                          cmpCovarConf,
                          cmpConf));
      } else {
        hybridBlockChain_.push_back(std::make_unique<SaberParametricBlockChain>(*hybridGeom,
                          *dualResGeom,
                          cmpOuterVars,
                          fset4dXb,
                          fset4dFg,
                          fset4dCmpEns,
                          fsetDualResEns,
                          cmpCovarConf,
                          cmpConf));
      }
    }
  } else {
    // Non-hybrid covariance: single block chain
    if (saberCentralBlockParams.saberBlockName.value() == "Ensemble") {
      hybridBlockChain_.push_back(std::make_unique<SaberEnsembleBlockChain>(geom,
                          *dualResGeom,
                          outerVars,
                          fset4dXb,
                          fset4dFg,
                          fsetEns,
                          fsetDualResEns,
                          covarConf,
                          params.toConfiguration()));
    } else {
      hybridBlockChain_.push_back(std::make_unique<SaberParametricBlockChain>(geom,
                          *dualResGeom,
                          outerVars,
                          fset4dXb,
                          fset4dFg,
                          fsetEns,
                          fsetDualResEns,
                          covarConf,
                          params.toConfiguration()));
    }

    // Set weights
    hybridScalarWeightSqrt_.push_back(1.0);
    // File-base weight
    atlas::FieldSet fsetWeight;
    hybridFieldWeightSqrt_.push_back(fsetWeight);
  }

  oops::Log::trace() << "ErrorCovariance::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::linearize(const State_ &xb,
                                       const Geometry_ &geom) {
  oops::Log::trace() << "ErrorCovariance linearized." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
size_t ErrorCovariance<MODEL>::ctlVecSize() const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::ctlVecSize starting" << std::endl;

  size_t ctlVecSize = 0;
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Add component
    ctlVecSize += hybridBlockChain_[jj]->ctlVecSize();
  }
  oops::Log::trace() << "ErrorCovariance<MODEL>::ctlVecSize done" << std::endl;
  return ctlVecSize;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::randomize(Increment_ & dx3d) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // JEDI compatibility
  Increment4D_ dx(dx3d.geometry(), dx3d.variables(), {dx3d.validTime()});

  // Loop over components for the central block
  oops::FieldSet4D fset4dSum(dx.times(), dx.commTime(), dx.geometry().getComm());
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Randomize covariance
    oops::FieldSet4D fset4dCmp(dx.times(), dx.commTime(), dx.geometry().getComm());
    hybridBlockChain_[jj]->randomize(fset4dCmp);

    // Weight square-root multiplication
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[jj];
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[jj];
    }

    if (jj == 0) {
      // Initialize sum
      for (size_t jtime = 0; jtime < fset4dSum.size(); ++jtime) {
        fset4dSum[jtime].fieldSet() = util::copyFieldSet(fset4dCmp[jtime].fieldSet());
      }
    } else {
      // Add component
      fset4dSum += fset4dCmp;
    }
  }

  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4dSum);

  // ATLAS fieldset to Increment_
  for (size_t jtime = 0; jtime < dx.size(); ++jtime) {
    dx[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  // JEDI compatibility
  dx3d = dx[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::multiply(const Increment_ &dx3di, Increment_ &dx3do) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // JEDI compatibility
  Increment4D_ dxi(dx3di.geometry(), dx3di.variables(), {dx3di.validTime()});
  Increment4D_ dxo(dx3do.geometry(), dx3do.variables(), {dx3do.validTime()});
  dxi[0] = dx3di;

  // Copy input
  dxo = dxi;
  oops::FieldSet4D fset4dInit(dxi);

  // Apply outer blocks adjoint
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocksAD(fset4dInit);

  // Initialize sum to zero
  oops::FieldSet4D fset4dSum = oops::copyFieldSet4D(fset4dInit);
  fset4dSum.zero();

  // Loop over B components
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Create temporary FieldSet
    oops::FieldSet4D fset4dCmp = oops::copyFieldSet4D(fset4dInit);

    // Apply weight
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[jj];
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[jj];
    }

    // Apply covariance
    hybridBlockChain_[jj]->multiply(fset4dCmp);

    // Apply weight
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[jj];
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[jj];
    }

    // Add component
    fset4dSum += fset4dCmp;
  }

  // Apply outer blocks forward
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4dSum);

  // ATLAS fieldset to Increment_
  for (size_t jtime = 0; jtime < dxo.size(); ++jtime) {
    dxo[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  // JEDI compatibility
  dx3do = dxo[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::inverseMultiply(const Increment_ &dxi, Increment_ &dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");

  // Iterative inverse
  oops::IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

  oops::Log::trace() << "ErrorCovariance<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::multiplySqrt(const IncrCtlVec_ &dv,
                                          Increment_ &dx3d) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrt starting" << std::endl;

  // JEDI compatibility
  Increment4D_ dx(dx3d.geometry(), dx3d.variables(), {dx3d.validTime()});

  // Loop over components for the central block
  size_t offset = 0;
  oops::FieldSet4D fset4dSum(dx.times(), dx.commTime(), dx.geometry().getComm());
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Apply covariance square-root
    oops::FieldSet4D fset4dCmp(dx.times(), dx.commTime(), dx.geometry().getComm());
    hybridBlockChain_[jj]->multiplySqrt(dv.modCtlVec().genCtlVec().data(), fset4dCmp, offset);
    offset += hybridBlockChain_[jj]->ctlVecSize();

    // Apply weight
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[jj];
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[jj];
    }

    if (jj == 0) {
      // Initialize sum
      for (size_t jtime = 0; jtime < fset4dSum.size(); ++jtime) {
        fset4dSum[jtime].fieldSet() = util::copyFieldSet(fset4dCmp[jtime].fieldSet());
      }
    } else {
      // Add component
      fset4dSum += fset4dCmp;
    }
  }

  // Apply outer blocks forward
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4dSum);

  // ATLAS fieldset to Increment_
  for (size_t jtime = 0; jtime < dx.size(); ++jtime) {
    dx[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  // JEDI compatibility
  dx3d = dx[0];

  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ErrorCovariance<MODEL>::multiplySqrtTrans(const Increment_ &dx3d,
                                               IncrCtlVec_ &dv) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrtTrans starting" << std::endl;

  // JEDI compatibility
  Increment4D_ dx(dx3d.geometry(), dx3d.variables(), {dx3d.validTime()});
  dx[0] = dx3d;

  // Create input FieldSet
  oops::FieldSet4D fset4dInit(dx);

  // Apply outer blocks adjoint
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocksAD(fset4dInit);

  // Loop over B components
  size_t offset = 0;
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Create temporary FieldSet
    oops::FieldSet4D fset4dCmp = oops::copyFieldSet4D(fset4dInit);

    // Apply weight
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[jj];
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[jj];
    }

    // Apply covariance square-root adjoint
    hybridBlockChain_[jj]->multiplySqrtTrans(fset4dCmp, dv.modCtlVec().genCtlVec().data(), offset);
    offset += hybridBlockChain_[jj]->ctlVecSize();
  }

  oops::Log::trace() << "ErrorCovariance<MODEL>::multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
