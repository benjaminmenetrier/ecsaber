/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/interface/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/interface/Increment.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/oops/ECUtilities.h"

namespace oops {
  class FieldSet3D;
}

namespace saber {

// -----------------------------------------------------------------------------

class SaberOuterBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberOuterBlockBase(const SaberBlockParametersBase & params)
    : blockName_(params.saberBlockName), skipInverse_(params.skipInverse),
      filterMode_(params.filterMode) {}
  virtual ~SaberOuterBlockBase() {}

  // Accessor

  // To inner Geometry data
  virtual const oops::GeometryData & innerGeometryData() const = 0;

  // To inner variables
  virtual const oops::patch::Variables & innerVars() const = 0;

  // Application methods

  // Block multiplication
  virtual void multiply(atlas::FieldSet &) const = 0;

  // Block multiplication adjoint
  virtual void multiplyAD(atlas::FieldSet &) const = 0;

  // Block left inverse multiplication
  virtual void leftInverseMultiply(atlas::FieldSet &) const

    {throw eckit::NotImplemented("leftInverseMultiply not implemented yet for the block "
      + blockName_, Here());}

  // Setup / calibration methods

  // Read block data
  virtual void read()
    {throw eckit::NotImplemented("read not implemented yet for the block " + this->blockName());}

  // Read model fields
  virtual std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToRead()
    {return {};}

  // Direct calibration
  virtual void directCalibration(const std::vector<atlas::FieldSet> &)
    {throw eckit::NotImplemented("directCalibration not implemented yet for the block "
      + this->blockName(), Here());}

  // Iterative calibration
  virtual void iterativeCalibrationInit()
    {throw eckit::NotImplemented("iterativeCalibrationInit not implemented yet for the block "
      + this->blockName(), Here());}
  virtual void iterativeCalibrationUpdate(const atlas::FieldSet &)
    {throw eckit::NotImplemented("iterativeCalibrationUpdate not implemented yet for the block "
      + this->blockName(), Here());}
  virtual void iterativeCalibrationFinal()
    {throw eckit::NotImplemented("iterativeCalibrationUpdate not implemented yet for the block "
      + this->blockName(), Here());}

  // Dual resolution setup
  virtual void dualResolutionSetup(const oops::GeometryData &)
    {throw eckit::NotImplemented("dualResolutionSetup not implemented yet for the block "
      + this->blockName(), Here());}

  // Write block data
  virtual void write() const {}

  // Write model fields
  virtual std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite() const
     {return {};}

  // Generate inner FieldSet (for the inverse test)
  virtual atlas::FieldSet generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                                const oops::patch::Variables & innerVars) const
    {return util::createRandomFieldSet(innerGeometryData.comm(),
                                       innerGeometryData.functionSpace(),
                                       innerVars);}

  // Generate outer FieldSet (for the inverse test)
  virtual atlas::FieldSet generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                                const oops::patch::Variables & outerVars) const
    {return util::createRandomFieldSet(outerGeometryData.comm(),
                                       outerGeometryData.functionSpace(),
                                       outerVars);}

  // Compare FieldSets (for the inverse test)
  virtual bool compareFieldSets(const atlas::FieldSet & fset1,
                                const atlas::FieldSet & fset2,
                                const double & tol) const
    {return util::compareFieldSets(fset1, fset2, tol);}

  // Non-virtual methods

  // Return block name
  std::string blockName() const {return blockName_;}

  // Return flag to skip inverse application
  bool filterMode() const {return filterMode_;}

  // Return flag to skip inverse application
  bool skipInverse() const {return skipInverse_;}

  // Read model fields
  template <typename MODEL>
  void read(const oops::Geometry<MODEL> &,
            const oops::patch::Variables &,
            const util::DateTime &);

  // Write model fields
  template <typename MODEL>
  void write(const oops::Geometry<MODEL> &,
             const oops::patch::Variables &,
             const util::DateTime &) const;

  // Adjoint test
  void adjointTest(const oops::GeometryData &,
                   const oops::patch::Variables &,
                   const oops::GeometryData &,
                   const oops::patch::Variables &,
                   const double &) const;

  // Left-inverse test
  void inverseTest(const oops::GeometryData &,
                   const oops::patch::Variables &,
                   const oops::GeometryData &,
                   const oops::patch::Variables &,
                   const oops::patch::Variables &,
                   const oops::patch::Variables &,
                   const double &,
                   const double &) const;

 private:
  std::string blockName_;
  bool skipInverse_;
  bool filterMode_;
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

class SaberOuterBlockFactory;

// -----------------------------------------------------------------------------

class SaberOuterBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberOuterBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberOuterBlockFactory>
    saberOuterBlockParameters{"saber block name", this};
};

// -----------------------------------------------------------------------------

class SaberOuterBlockFactory {
 public:
  static std::unique_ptr<SaberOuterBlockBase> create(const oops::GeometryData &,
                                                     const oops::patch::Variables &,
                                                     const eckit::Configuration &,
                                                     const SaberBlockParametersBase &,
                                                     const oops::FieldSet3D &,
                                                     const oops::FieldSet3D &);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberOuterBlockFactory() = default;

 protected:
  explicit SaberOuterBlockFactory(const std::string &name);

 private:
  virtual std::unique_ptr<SaberOuterBlockBase> make(const oops::GeometryData &,
                                                    const oops::patch::Variables &,
                                                    const eckit::Configuration &,
                                                    const SaberBlockParametersBase &,
                                                    const oops::FieldSet3D &,
                                                    const oops::FieldSet3D &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberOuterBlockFactory * > & getMakers() {
    static std::map < std::string, SaberOuterBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberOuterBlockMaker : public SaberOuterBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<SaberOuterBlockBase> make(const oops::GeometryData & outerGeometryData,
                                            const oops::patch::Variables & outerVars,
                                            const eckit::Configuration & covarConf,
                                            const SaberBlockParametersBase & params,
                                            const oops::FieldSet3D & xb,
                                            const oops::FieldSet3D & fg) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::make_unique<T>(outerGeometryData, outerVars,
                               covarConf, stronglyTypedParams, xb, fg);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberOuterBlockMaker(const std::string & name) : SaberOuterBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterBlockBase::read(const oops::Geometry<MODEL> & geom,
                               const oops::patch::Variables & vars,
                               const util::DateTime & date) {
  oops::Log::trace() << "SaberOuterBlockBase::read starting" << std::endl;

  // Get vector of configuration/FieldSet pairs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> inputs
    = this->fieldsToRead();

  // Read fieldsets as increments
  for (auto & input : inputs) {
    // Create increment
    oops::Increment<MODEL> dx(geom, templatedVars<MODEL>(vars), date);
    dx.read(input.first);
    oops::Log::test() << "Norm of input parameter " << input.second.name()
                      << ": " << dx.norm() << std::endl;
    util::copyFieldSet(dx.increment().fieldSet(), input.second);
  }

  oops::Log::trace() << "SaberOuterBlockBase::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterBlockBase::write(const oops::Geometry<MODEL> & geom,
                                const oops::patch::Variables & vars,
                                const util::DateTime & date) const {
  oops::Log::trace() << "SaberOuterBlockBase::write starting" << std::endl;

  // Get vector of configuration/FieldSet pairs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> outputs
    = this->fieldsToWrite();

  // Create increment
  oops::Increment<MODEL> dx(geom, templatedVars<MODEL>(vars), date);

  // Loop and write
  for (const auto & output : outputs) {
    dx.increment().fieldSet() = util::copyFieldSet(output.second);
    dx.increment().synchronizeFields();
    oops::Log::test() << "Norm of output parameter " << output.second.name()
                      << ": " << dx.norm() << std::endl;
    dx.write(output.first);
  }

  oops::Log::trace() << "SaberOuterBlockBase::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber