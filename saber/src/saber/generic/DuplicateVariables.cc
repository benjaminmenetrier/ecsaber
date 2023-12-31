/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/DuplicateVariables.h"

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

namespace {

oops::patch::Variables createActiveVars(const std::vector<VariableGroupParameters> & gps,
                                 const oops::patch::Variables & outerVars)  {
  oops::patch::Variables activeVars;
  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    oops::patch::Variables v = gp.groupComponents.value();
    activeVars += v;
    for (const std::string & s : v.variables()) {
      activeVars.addMetaData(s, "levels", outerVars.getLevels(s));
    }
    activeVars.push_back(key);
    activeVars.addMetaData(key, "levels", outerVars.getLevels(v[0]));
  }

  return activeVars;
}


oops::patch::Variables createInnerVars(const oops::patch::Variables & outerVars,
                                const oops::patch::Variables & activeVars,
                                const std::vector<VariableGroupParameters> & gps) {
  oops::patch::Variables innerVars;
  for (const std::string & var : outerVars.variables()) {
    if (!activeVars.has(var)) {
      innerVars.push_back(var);
      innerVars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
  }

  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    innerVars.push_back(key);
    innerVars.addMetaData(key, "levels", activeVars.getLevels(key));
  }

  return innerVars;
}

// copy groupVars fields (assuming all activeFields are allocated)
void copyFields(const std::vector<VariableGroupParameters> & gps,
                const atlas::FieldSet & fsetIn,
                atlas::FieldSet & fsetOut) {
  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    oops::patch::Variables v = gp.groupComponents.value();
    auto otherView = atlas::array::make_view<double, 2>(fsetIn[key]);
    for (const std::string & component : v.variables()) {
      auto view = atlas::array::make_view<double, 2>(fsetOut[component]);
      view.assign(otherView);
    }
  }
}

// gather
void gatherFields(const std::vector<VariableGroupParameters> & gps,
                  atlas::FieldSet & fsetIn,
                  atlas::FieldSet & fsetOut) {
  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    oops::patch::Variables v = gp.groupComponents.value();
    auto otherView = atlas::array::make_view<double, 2>(fsetOut[key]);
    for (const std::string & component : v.variables()) {
      auto view = atlas::array::make_view<double, 2>(fsetIn[component]);
      for (atlas::idx_t jn = 0; jn < fsetOut[key].shape(0); ++jn) {
        for (atlas::idx_t jl = 0; jl < fsetOut[key].shape(1); ++jl) {
          otherView(jn, jl) += view(jn, jl);
          view(jn, jl) = 0.0;
        }
      }
    }
  }
}

}  //  namespace


// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DuplicateVariables>
    makerDuplicateVariables_("duplicate variables");

// -----------------------------------------------------------------------------

DuplicateVariables::DuplicateVariables(const oops::GeometryData & outerGeometryData,
                                       const oops::patch::Variables & outerVars,
                                       const eckit::Configuration & covarConfig,
                                       const Parameters_ & params,
                                       const oops::FieldSet3D & xb,
                                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    groups_(params.variableGroupParameters.value()),
    outerVars_(outerVars),
    activeVars_(createActiveVars(groups_, outerVars_)),
    innerVars_(createInnerVars(outerVars_, activeVars_, groups_)),
    innerGeometryData_(outerGeometryData)
{
  oops::Log::trace() << classname() << "::DuplicateVariables starting" << std::endl;
  oops::Log::trace() << classname() << "::DuplicateVariables done" << std::endl;
}

// -----------------------------------------------------------------------------
void DuplicateVariables::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // The aim is to copy groupVars_ fields into component fields.
  // Some manipulation of fieldsets is needed for this.
  atlas::FieldSet fsetOut = atlas::FieldSet();

  for (const VariableGroupParameters & gp : groups_) {
    oops::patch::Variables v = gp.groupComponents.value();
    for (const std::string& component : v.variables()) {
      const size_t nlev = activeVars_.getLevels(component);
      atlas::Field field =
        innerGeometryData_.functionSpace()->createField<double>(
        atlas::option::name(component) | atlas::option::levels(nlev));
      atlas::array::make_view<double, 2>(field).assign(0.0);
      if (innerGeometryData_.functionSpace()->type() != "Spectral") {
        field.haloExchange();
      }
      fsetOut.add(field);
    }
  }

  // copy group fields into component fields
  copyFields(groups_, fset.fieldSet(), fsetOut);

  // keep passive vars
  for (auto & fld : fset) {
    if (!activeVars_.has(fld.name())) {
      fsetOut.add(fld);
    }
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DuplicateVariables::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // The aim is to do the adjoint of the copy groupVars_ fields into component fields.
  // This involves summing at the component fields together in
  // Some manipulation of fieldsets is needed for this.
  atlas::FieldSet fsetOut = atlas::FieldSet();

  // allocate group vars.
  for (const VariableGroupParameters & gp : groups_) {
    std::string key = gp.groupVariableName.value();
    const size_t nlev = activeVars_.getLevels(key);
    atlas::Field field =
      innerGeometryData_.functionSpace()->createField<double>(
        atlas::option::name(key) | atlas::option::levels(nlev));
    atlas::array::make_view<double, 2>(field).assign(0.0);
    if (innerGeometryData_.functionSpace()->type() != "Spectral") {
      field.haloExchange();
    }
    fsetOut.add(field);
  }

  // sum component fields into group fields
  gatherFields(groups_, fset.fieldSet(), fsetOut);

  // keep passive vars
  for (auto & fld : fset) {
    if (!activeVars_.has(fld.name())) {
      fsetOut.add(fld);
    }
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiplyAD done = " << std::endl;
}

// -----------------------------------------------------------------------------

void DuplicateVariables::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
