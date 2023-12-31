/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/covariance/Covariance.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/gsi/covariance/Covariance.interface.h"
#include "saber/gsi/grid/Grid.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberCentralBlockMaker<Covariance> makerCovariance_("gsi covariance");

// -------------------------------------------------------------------------------------------------

Covariance::Covariance(const oops::GeometryData & geometryData,
                       const oops::patch::Variables & centralVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime()),
    params_(params), variables_(params.activeVars.value().get_value_or(centralVars).variables()),
    gsiGridFuncSpace_(geometryData.functionSpace()), comm_(&geometryData.comm()),
    xb_(xb), fg_(fg), validTimeOfXbFg_(xb.validTime())
{
  oops::Log::trace() << classname() << "::Covariance starting" << std::endl;
  oops::Log::trace() << classname() << "::Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Covariance::~Covariance() {
  oops::Log::trace() << classname() << "::~Covariance starting" << std::endl;
  util::Timer timer(classname(), "~Covariance");
  gsi_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // Ignore incoming fields and create new ones based on the block function space
  // ----------------------------------------------------------------------------
  atlas::FieldSet newFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding fields on gsi grid
  for (const auto & sabField : fset) {
      // Ensure that the field name is in the variables list
      if (std::find(variables_.begin(), variables_.end(), sabField.name()) == variables_.end()) {
        throw eckit::Exception("Field " + sabField.name() + " not found in the " + classname() +
          " variables.", Here());
      }

      // Create the gsi grid field and add to Fieldset
      newFields.add(gsiGridFuncSpace_.createField<double>(atlas::option::name(sabField.name())
        | atlas::option::levels(sabField.levels())));
  }

  // Replace whatever fields are coming in with the gsi grid fields
  fset.fieldSet() = newFields;

  // Call implementation
  gsi_covariance_randomize_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  gsi_covariance_multiply_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Covariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Create covariance module
  gsi_covariance_create_f90(keySelf_, *comm_, params_.readParams.value()->toConfiguration(),
     xb_.get(), fg_.get(), validTimeOfXbFg_);
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
