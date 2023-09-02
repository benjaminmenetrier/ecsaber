/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/generic/GenericCtlVec.h"
#include "oops/generic/LocalizationBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template<typename MODEL>
class Localization : public oops::LocalizationBase<MODEL> {
  using GenericCtlVec_ = oops::GenericCtlVec;
  using Geometry_ = oops::Geometry<MODEL>;
  using Increment_ = oops::Increment<MODEL>;
  using Model_ = oops::Model<MODEL>;
  using State_ = oops::State<MODEL>;

 public:
  Localization(const Geometry_ &,
               const oops::Variables &,
               const eckit::Configuration &);
  ~Localization();

  void multiply(Increment_ &) const override;

  // Square-root formulation
  size_t ctlVecSize() const {return loc_->ctlVecSize();}
  void multiplySqrt(const GenericCtlVec_ & dv, Increment_ & dx) const;
  void multiplySqrtTrans(const Increment_ & dx, GenericCtlVec_ & dv) const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<SaberParametricBlockChain> loc_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & geom,
                                  const oops::Variables & incVarsNoMeta,
                                  const eckit::Configuration & conf)
  : loc_()
{
  oops::Log::trace() << "Localization::Localization starting" << std::endl;

  // Create dummy time
  util::DateTime dummyTime(1977, 5, 25, 0, 0, 0);

  // Initialize
  const std::vector<std::size_t> vlevs = geom.variableSizes(incVarsNoMeta);
  oops::Variables incVars(incVarsNoMeta);
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    incVars.addMetaData(incVarsNoMeta.variables()[i], "levels", vlevs[i]);
  }

  // Create dummy xb and fg
  atlas::FieldSet fset = util::createRandomFieldSet(geom.generic().comm(),
                                                    geom.functionSpace(),
                                                    incVars);
  oops::FieldSet3D fset3d(fset, dummyTime, eckit::mpi::comm());
  oops::FieldSet4D fset4dXb(fset3d);
  oops::FieldSet4D fset4dFg(fset3d);

  std::vector<atlas::FieldSet> fsetEns;
  // TODO(AS): revisit what configuration needs to be passed to SaberParametricBlockChain.
  eckit::LocalConfiguration covarConf;
  eckit::LocalConfiguration ensembleConf;
  ensembleConf.set("ensemble size", 0);
  covarConf.set("ensemble configuration", ensembleConf);
  covarConf.set("adjoint test", conf.getBool("adjoint test", false));
  covarConf.set("adjoint tolerance", conf.getDouble("adjoint tolerance", 1.0e-12));
  covarConf.set("inverse test", conf.getBool("inverse test", false));
  covarConf.set("inverse tolerance", conf.getDouble("inverse tolerance", 1.0e-12));
  covarConf.set("square-root test", conf.getBool("square-root test", false));
  covarConf.set("square-root tolerance", conf.getDouble("square-root tolerance", 1.0e-12));
  covarConf.set("iterative ensemble loading", false);
  // Initialize localization blockchain
  loc_ = std::make_unique<SaberParametricBlockChain>(geom, geom,
              incVars, fset4dXb, fset4dFg,
              fsetEns, fsetEns, covarConf, conf);

  oops::Log::trace() << "Localization:Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  oops::Log::trace() << "Localization:~Localization destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  oops::Log::trace() << "Localization:multiply starting" << std::endl;

  // SABER block chain multiplication
  oops::FieldSet4D fset4d({dx.fieldSet(), dx.validTime(), dx.geometry().getComm()});
  loc_->multiply(fset4d);

  // ATLAS fieldset to Increment_
  dx.synchronizeFields();

  oops::Log::trace() << "Localization:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiplySqrt(const GenericCtlVec_ & dv,
                                       Increment_ & dx) const {
  oops::Log::trace() << "Localization:multiplySqrt starting" << std::endl;

  // SABER block chain square-root
  oops::FieldSet3D fset3d(dx.validTime(), dx.geometry().getComm());
  oops::FieldSet4D fset4d(fset3d);
  loc_->multiplySqrt(dv.data(), fset4d, 0);

  // ATLAS fieldset to Increment_
  dx.fromFieldSet(fset4d[0].fieldSet());

  oops::Log::trace() << "Localization:multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiplySqrtTrans(const Increment_ & dx,
                                            GenericCtlVec_ & dv) const {
  oops::Log::trace() << "Localization:multiplySqrtTrans starting" << std::endl;

  // SABER block chain square-root adjoint
  oops::FieldSet4D fset4d({dx.fieldSet(), dx.validTime(), dx.geometry().getComm()});
  loc_->multiplySqrtTrans(fset4d, dv.data(), 0);

  oops::Log::trace() << "Localization:multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  os << "Localization:print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace saber
