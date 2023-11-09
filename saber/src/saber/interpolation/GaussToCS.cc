/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "oops/util/Logger.h"

#include "saber/interpolation/GaussToCS.h"

using atlas::grid::detail::partitioner::TransPartitioner;
using atlas::grid::detail::partitioner::MatchingMeshPartitionerCubedSphere;

namespace saber {
namespace interpolation {

namespace {

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid) {
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(1));
}

// -----------------------------------------------------------------------------

auto createPointCloud(const atlas::Grid& grid,
                      const atlas::grid::Partitioner& partitioner) {
    const auto distribution = atlas::grid::Distribution(grid, partitioner);

    auto lonLats = std::vector<atlas::PointXY>{};
    auto idx = atlas::gidx_t{0};
    for (const auto& lonLat : grid.lonlat()) {
        if (distribution.partition(idx++) == atlas::mpi::rank()) {
            lonLats.emplace_back(lonLat.data());
        }
    }
    return std::make_unique<atlas::functionspace::PointCloud>(
                atlas::functionspace::PointCloud(lonLats));
}

// -----------------------------------------------------------------------------

auto createInverseInterpolation(const bool initializeInverseInterpolation,
                                const bool hasSinglePE,
                                const atlas::functionspace::NodeColumns & csFunctionSpace,
                                const atlas::Grid & gaussGrid,
                                const atlas::grid::Partitioner & gaussPartitioner) {
  CS2Gauss inverseInterpolation;

  if (!initializeInverseInterpolation || hasSinglePE) {
    return inverseInterpolation;
  }

  const auto config = atlas::util::Config("type", "cubedsphere");
  const atlas::grid::MatchingMeshPartitioner csMatchingPartitioner(csFunctionSpace.mesh(),
                                                                   config);
  try {
    inverseInterpolation.matchingPtcldFspace = createPointCloud(gaussGrid, csMatchingPartitioner);
  } catch (const eckit::Exception &) {
    oops::Log::error()
            << "ERROR  : SABER block \"gauss to cubed-sphere dual\" inverse cannot run"
            << std::endl
            << "ERROR  : with 2 or 3 MPI tasks from an " << gaussGrid.name()
            << " grid to a " << csFunctionSpace.mesh().grid().name() << " grid."
            << std::endl
            << "ERROR  : Try one of the three possible solutions:" << std::endl
            << "ERROR  :   1. Add yaml key `initialize inverse interpolator: false` to the block"
            << std::endl
            << "ERROR  :   2. Use 1, 4 or more MPI tasks" << std::endl
            << "ERROR  :   3. Use a Gaussian grid with less points" << std::endl;
    throw(eckit::FunctionalityNotSupported(
                "Inverse not implemented on 2 and 3 PEs with these grids.", Here()));
  }

  inverseInterpolation.targetPtcldFspace = createPointCloud(gaussGrid, gaussPartitioner);

  const auto interpConfig =  atlas::util::Config("type", "cubedsphere-bilinear") |
                             atlas::util::Config("halo_exchange", false);

  inverseInterpolation.interpolation = atlas::Interpolation(
              interpConfig, csFunctionSpace,
              *inverseInterpolation.matchingPtcldFspace);

  inverseInterpolation.redistribution = atlas::Redistribution(
              *inverseInterpolation.matchingPtcldFspace,
              *inverseInterpolation.targetPtcldFspace);

  return inverseInterpolation;
}

// -----------------------------------------------------------------------------

/* Direct interpolation from NodeColumn cubed-sphere FunctionSpace to
 * StructuredColumns is not possible on multiple PEs.
 * Here, the interpolation is done following this route:
 * CubedSphere FunctionSpace
 * -> matching PointCloud FunctionSpace using CS partitioner
 * -> target PointCloud FunctionSpace using Gauss grid partitioner
 * -> gauss StructuredColumns FunctionSpace.
 */
void inverseInterpolateMultiplePEs(
        const oops::patch::Variables & variables,
        const CS2Gauss & inverseInterpolation,
        const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
        const atlas::FieldSet & srcFieldSet,
        atlas::FieldSet & newFieldSet) {
  // Interpolate from source to matching PointCloud
  atlas::FieldSet matchingPtcldFset;
  for (const auto & fieldname : variables.variables()) {
    auto matchingPtcldField = inverseInterpolation.matchingPtcldFspace->createField<double>(
                atlas::option::name(fieldname) |
                atlas::option::levels(srcFieldSet[fieldname].levels()));
    atlas::array::make_view<double, 2>(matchingPtcldField).assign(0.0);
    matchingPtcldFset.add(matchingPtcldField);
  }
  inverseInterpolation.interpolation.execute(srcFieldSet, matchingPtcldFset);

  // Redistribute from matching PointCloud to target PointCloud
  atlas::FieldSet targetPtcldFset;
  for (const auto & fieldname : variables.variables()) {
    auto targetPtcldField = inverseInterpolation.targetPtcldFspace->createField<double>(
          atlas::option::name(fieldname) |
          atlas::option::levels(srcFieldSet[fieldname].levels()));
    targetPtcldFset.add(targetPtcldField);
  }
  inverseInterpolation.redistribution.execute(matchingPtcldFset, targetPtcldFset);

  // Copy from target PointCloud to gauss StructuredColumns
  for (const auto & fieldname : variables.variables()) {
    atlas::Field gaussField = gaussFunctionSpace.createField<double>(
                atlas::option::name(fieldname) |
                atlas::option::levels(srcFieldSet[fieldname].levels()));
    atlas::array::make_view<double, 2>(gaussField).assign(
        atlas::array::make_view<const double, 2>(targetPtcldFset[fieldname]));
    gaussField.haloExchange();
    newFieldSet.add(gaussField);
  }
}

// -----------------------------------------------------------------------------

void inverseInterpolateSinglePE(
        const oops::patch::Variables & variables,
        const atlas::functionspace::NodeColumns & CSFunctionSpace,
        const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
        const atlas::Grid & gaussGrid,
        const atlas::FieldSet & srcFieldSet,
        atlas::FieldSet & newFieldSet) {
  const auto targetMesh = atlas::MeshGenerator("structured").generate(gaussGrid);
  const auto hybridFunctionSpace = atlas::functionspace::NodeColumns(targetMesh);
  const auto scheme = atlas::util::Config("type", "cubedsphere-bilinear") |
                      atlas::util::Config("adjoint", "false");
  const auto interp = atlas::Interpolation(scheme, CSFunctionSpace, hybridFunctionSpace);

  atlas::FieldSet hybridFieldSet;
  for (const auto & fieldname : variables.variables()) {
    atlas::Field hybridField = hybridFunctionSpace.createField<double>(
                atlas::option::name(fieldname) |
                atlas::option::levels(srcFieldSet[fieldname].levels()));
    hybridFieldSet.add(hybridField);
  }

  interp.execute(srcFieldSet, hybridFieldSet);

  // Copy into StructuredColumns
  for (const auto & fieldname : variables.variables()) {
    atlas::Field gaussField = gaussFunctionSpace.createField<double>(
                atlas::option::name(fieldname) |
                atlas::option::levels(srcFieldSet[fieldname].levels()));
    atlas::array::make_view<double, 2>(gaussField).assign(
        atlas::array::make_view<const double, 2>(hybridFieldSet[fieldname]));
    gaussField.haloExchange();
    newFieldSet.add(gaussField);
  }
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussToCS> makerGaussToCS_("gauss to cubed-sphere-dual");

// -----------------------------------------------------------------------------
// Note that this is slower than this needs to be
// as we are need to create 2 grid objects (very slow)
// In the future it might make sense to include an atlas grid (if available) from
// the model in outerGeometryData.
GaussToCS::GaussToCS(const oops::GeometryData & outerGeometryData,
                     const oops::patch::Variables & outerVars,
                     const eckit::Configuration & covarConf,
                     const Parameters_ & params,
                     const oops::FieldSet3D & xb,
                     const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params),
    innerVars_(outerVars),
    activeVars_(params.activeVariables.value().get_value_or(innerVars_)),
    CSFunctionSpace_(outerGeometryData.functionSpace()),
    gaussGrid_(params.gaussGridUid.value()),
    gaussFunctionSpace_(createGaussFunctionSpace(gaussGrid_)),
    gaussPartitioner_(new TransPartitioner()),
    csgrid_(CSFunctionSpace_.mesh().grid()),
    interp_(gaussPartitioner_, gaussFunctionSpace_, csgrid_, CSFunctionSpace_),
    inverseInterpolation_(createInverseInterpolation(
                              params.initializeInverseInterpolation.value(),
                              outerGeometryData.comm().size() == 1,
                              CSFunctionSpace_, gaussGrid_,
                              gaussPartitioner_)),
    innerGeometryData_(gaussFunctionSpace_, outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(),
                       outerGeometryData.comm())

{
  oops::Log::trace() << classname() << "::GaussToCS starting" << std::endl;
  oops::Log::trace() << classname() << "::GaussToCS done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet gaussFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       gaussFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // On input: fieldset on gaussian mesh

  // Create fieldset on cubed-sphere mesh.

  atlas::FieldSet csFieldSet;
  for (const auto & fieldname : activeVars_.variables()) {
    atlas::Field csField =
      CSFunctionSpace_.createField<double>(
          atlas::option::name(fieldname) |
          atlas::option::levels(gaussFieldSet[fieldname].levels()) |
          atlas::option::halo(1));
    csFieldSet.add(csField);
  }

  // Interpolate to cubed sphere
  interp_.execute(gaussFieldSet, csFieldSet);

  for (const auto & fieldname : activeVars_.variables()) {
    newFields.add(csFieldSet[fieldname]);
  }

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiply done"
                     << fieldSet.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // On input: fieldset on gaussian grid
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet csFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : fieldSet.field_names()) {
    if (activeVars_.has(fieldname)) {
      csFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // Create gauss fieldset
  atlas::FieldSet gaussFieldSet;
  for (const auto & fieldname : activeVars_.variables()) {
    atlas::Field gaussField =
      gaussFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
            atlas::option::levels(csFieldSet[fieldname].levels()) |
            atlas::option::halo(1));
    gaussFieldSet.add(gaussField);
  }

  // Adjoint of interpolation from gauss to dual cubed sphere
  interp_.executeAdjoint(gaussFieldSet, csFieldSet);

  for (const auto & fieldname : activeVars_.variables()) {
    newFields.add(gaussFieldSet[fieldname]);
  }

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::leftInverseMultiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  atlas::FieldSet newFieldSet = atlas::FieldSet();
  atlas::FieldSet srcFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       srcFieldSet.add(fieldSet[fieldname]);
     } else {
       newFieldSet.add(fieldSet[fieldname]);
     }
  }

  if (innerGeometryData_.comm().size() >= 2) {
    inverseInterpolateMultiplePEs(activeVars_, inverseInterpolation_,
                                 gaussFunctionSpace_,
                                 srcFieldSet, newFieldSet);
  } else {
    // Approach above uses the cubed-sphere matching partitioner, which fails
    // to partition on 1-3 PEs for some configurations (e.g. from a CS-LFR-12
    // to F15). As a temporary and partial fix, we use a bespoke interpolation
    // when running on one MPI task.
    inverseInterpolateSinglePE(activeVars_,
                               CSFunctionSpace_,
                               gaussFunctionSpace_,
                               gaussGrid_,
                               srcFieldSet, newFieldSet);
  }

  fieldSet = newFieldSet;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
