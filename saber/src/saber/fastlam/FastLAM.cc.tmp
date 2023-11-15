/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/FastLAM.h"

#include <math.h>
#include <netcdf.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<FastLAM> makerFastLAM_("FastLAM");

// -----------------------------------------------------------------------------

FastLAM::FastLAM(const oops::GeometryData & gdata,
                 const oops::patch::Variables & activeVars,
                 const eckit::Configuration & covarConf,
                 const Parameters_ & params,
                 const oops::FieldSet3D & xb,
                 const oops::FieldSet3D & fg) :
    SaberCentralBlockBase(params),
    gdata_(gdata),
    comm_(gdata_.comm()),
    activeVars_(activeVars),
    params_(params.calibration.value() != boost::none ? *params.calibration.value()
      : *params.read.value())
{
  oops::Log::trace() << classname() << "::FastLAM starting" << std::endl;

  // Check geometry fields
  ASSERT(gdata_.fieldSet().has_field("index_i"));
  ASSERT(gdata_.fieldSet().has_field("index_j"));
  ASSERT(gdata_.fieldSet().has_field("area"));

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Get grid size
  nx0_ = 0;
  ny0_ = 0;
  atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
  atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
  nodes0_ = fieldIndexI0.shape(0);
  auto viewIndexI0 = atlas::array::make_view<int, 1>(fieldIndexI0);
  auto viewIndexJ0 = atlas::array::make_view<int, 1>(fieldIndexJ0);
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    if (viewGhost0(jnode0) == 0) {
      nx0_ = std::max(nx0_, static_cast<size_t>(viewIndexI0(jnode0)));
      ny0_ = std::max(ny0_, static_cast<size_t>(viewIndexJ0(jnode0)));
    }
  }
  comm_.allReduceInPlace(nx0_, eckit::mpi::max());
  comm_.allReduceInPlace(ny0_, eckit::mpi::max());
  oops::Log::info() << "Regional grid size: " << nx0_ << "x" << ny0_ << std::endl;

  // Define 2d active variables
  active2dVars_ = oops::patch::Variables();
  for (const auto & var : activeVars_.variables()) {
    if (activeVars_.getLevels(var) == 1) {
      active2dVars_.push_back(var);
    }
  }

  // Create groups
  if (params_.groups.value() == boost::none) {
    // No group specified, each variable is its own group
    for (const auto & var : activeVars_.variables()) {
      // Define group properties
      const size_t groupNz0 = activeVars_.getLevels(var);
      const std::vector<std::string> variables = {var};

      // Add group
      groups_.push_back(std::make_tuple(var,
                                        groupNz0,
                                        var,
                                        variables));
    }
  } else {
    // Copy groups
    for (const auto & group : *params_.groups.value()) {
      // Define group properties
      const std::string groupName = group.name.value();
      size_t groupNz0 = 1;
      for (const auto & var : group.variables.value()) {
        if (!active2dVars_.has(var)) {
          if (groupNz0 == 1) {
            // Assign number of levels
            groupNz0 = static_cast<size_t>(activeVars_.getLevels(var));
          } else {
            // Check number of levels
            ASSERT(static_cast<int>(groupNz0) == activeVars_.getLevels(var));
          }
        }
      }

      // Add group
      groups_.push_back(std::make_tuple(groupName,
                                        groupNz0,
                                        group.varInModelFile.value().get_value_or(groupName),
                                        group.variables.value()));
    }
  }

  oops::Log::trace() << classname() << "::FastLAM done" << std::endl;
}

// -----------------------------------------------------------------------------

FastLAM::~FastLAM() {
  oops::Log::trace() << classname() << "::~FastLAM starting" << std::endl;
  oops::Log::trace() << classname() << "::~FastLAM done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("randomize method not implemented yet", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Save input FieldSet
  atlas::FieldSet fsetIn = util::copyFieldSet(fset);
  util::zeroFieldSet(fset);

  // Loop over bins
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    // Copy input FieldSet
    atlas::FieldSet fsetBin = util::copyFieldSet(fsetIn);

    if (params_.strategy.value() == "univariate") {
      // Univariate strategy
      for (const auto & var : activeVars_.variables()) {
        // Group properties
        const std::string groupName = std::get<0>(groups_[getGroupIndex(var)]);

        // Variable properties
        const size_t varNz0 = activeVars_.getLevels(var);
        const size_t k0Offset = getK0Offset(var);

        // Apply weight square-root and normalization
        atlas::Field fieldBin = fsetBin[var];
        auto viewBin = atlas::array::make_view<double, 2>(fieldBin);
        const atlas::Field fieldWgtSqrt = weight_[jBin][groupName];
        const auto viewWgtSqrt = atlas::array::make_view<double, 2>(fieldWgtSqrt);
        const atlas::Field fieldNorm = normalization_[jBin][groupName];
        const auto viewNorm = atlas::array::make_view<double, 2>(fieldNorm);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (viewGhost0(jnode0) == 0) {
            for (size_t k0 = 0; k0 < varNz0; ++k0) {
              viewBin(jnode0, k0) *= viewWgtSqrt(jnode0, k0Offset+k0)*viewNorm(jnode0, k0Offset+k0);
            }
          }
        }

        // Layer multiplication
        data_.at(groupName)[jBin].multiply(fieldBin);

        // Apply weight square-root and normalization
        viewBin = atlas::array::make_view<double, 2>(fieldBin);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (viewGhost0(jnode0) == 0) {
            for (size_t k0 = 0; k0 < varNz0; ++k0) {
              viewBin(jnode0, k0) *= viewWgtSqrt(jnode0, k0Offset+k0)*viewNorm(jnode0, k0Offset+k0);
            }
          }
        }
      }
    } else if (params_.strategy.value() == "duplicated") {
      // Duplicated strategy
      for (const auto & group : groups_) {
        // Group properties
        const std::string groupName = std::get<0>(group);
        const size_t groupNz0 = std::get<1>(group);

        // Create group field
        atlas::Field fieldGrp = gdata_.functionSpace().createField<double>(
            atlas::option::name(groupName) | atlas::option::levels(groupNz0));
        auto viewGrp = atlas::array::make_view<double, 2>(fieldGrp);
        viewGrp.assign(0.0);

        // Sum all variables of the group
        for (const auto & var : std::get<3>(group)) {
          // Variable properties
          const size_t varNz0 = activeVars_.getLevels(var);
          const size_t k0Offset = getK0Offset(var);

          // Add variable field
          const atlas::Field fieldBin = fsetBin[var];
          const auto viewBin = atlas::array::make_view<double, 2>(fieldBin);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (viewGhost0(jnode0) == 0) {
              for (size_t k0 = 0; k0 < varNz0; ++k0) {
                viewGrp(jnode0, k0Offset+k0) += viewBin(jnode0, k0);
              }
            }
          }
        }

        // Apply weight square-root and normalization
        const atlas::Field fieldWgtSqrt = weight_[jBin][groupName];
        const auto viewWgtSqrt = atlas::array::make_view<double, 2>(fieldWgtSqrt);
        const atlas::Field fieldNorm = normalization_[jBin][groupName];
        const auto viewNorm = atlas::array::make_view<double, 2>(fieldNorm);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (viewGhost0(jnode0) == 0) {
            for (size_t k0 = 0; k0 < groupNz0; ++k0) {
              viewGrp(jnode0, k0) *= viewWgtSqrt(jnode0, k0)*viewNorm(jnode0, k0);
            }
          }
        }

        // Layer multiplication
        data_.at(groupName)[jBin].multiply(fieldGrp);

        // Apply weight square-root and normalization
        viewGrp = atlas::array::make_view<double, 2>(fieldGrp);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (viewGhost0(jnode0) == 0) {
            for (size_t k0 = 0; k0 < groupNz0; ++k0) {
              viewGrp(jnode0, k0) *= viewWgtSqrt(jnode0, k0)*viewNorm(jnode0, k0);
            }
          }
        }

        // Copy result on all variables of the group
        for (const auto & var : std::get<3>(group)) {
          // Variable properties
          const size_t varNz0 = activeVars_.getLevels(var);
          const size_t k0Offset = getK0Offset(var);

          // Add variable field
          atlas::Field fieldBin = fsetBin[var];
          auto viewBin = atlas::array::make_view<double, 2>(fieldBin);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (viewGhost0(jnode0) == 0) {
              for (size_t k0 = 0; k0 < varNz0; ++k0) {
                viewBin(jnode0, k0) = viewGrp(jnode0, k0Offset+k0);
              }
            }
          }
        }
      }
    } else {
      // Wrong multivariate strategy
      throw eckit::UserError("wrong multivariate strategy: " + params_.strategy.value(), Here());
    }

    // Add component
    util::addFieldSets(fset, fsetBin);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> FastLAM::fieldsToRead() {
  oops::Log::trace() << classname() << "::fieldsToRead starting" << std::endl;

  if (params_.inputModelFilesConf.value() != boost::none) {
    for (const auto & conf : *params_.inputModelFilesConf.value()) {
      // Get file configuration
      eckit::LocalConfiguration file = getFileConf(comm_, conf);

      // Get parameter
      const std::string param = conf.getString("parameter");

      // Number of components
      size_t nWeight = 0;
      size_t nNormalization = 0;

      if (param == "weight") {
        // Specific case for weight
        nWeight = conf.getInt("number of components");
        for (size_t jBin = 0; jBin < nWeight; ++jBin) {
          // Create FieldSet
          atlas::FieldSet fset;
          fset.name() = param + "_" + std::to_string(jBin);

          // Update configuration
          eckit::LocalConfiguration fileCmp(file);
          util::seekAndReplace(fileCmp, "%component%", jBin, 2);

          // Add pair
          inputs_.push_back(std::make_pair(fileCmp, fset));
        }
      } else if (param == "normalization") {
        // Specific case for normalization
        nNormalization = conf.getInt("number of components");
        for (size_t jBin = 0; jBin < nNormalization; ++jBin) {
          // Create FieldSet
          atlas::FieldSet fset;
          fset.name() = param + "_" + std::to_string(jBin);

          // Update configuration
          eckit::LocalConfiguration fileCmp(file);
          util::seekAndReplace(fileCmp, "%component%", jBin, 2);

          // Add pair
          inputs_.push_back(std::make_pair(fileCmp, fset));
        }
      } else {
        // Create FieldSet
        atlas::FieldSet fset;
        fset.name() = param;

        // Add pair
        inputs_.push_back(std::make_pair(file, fset));
      }
    }
  }

  oops::Log::trace() << classname() << "::fieldsToRead done" << std::endl;
  return inputs_;
}

// -----------------------------------------------------------------------------

void FastLAM::directCalibration(const std::vector<atlas::FieldSet> &) {
  oops::Log::trace() << classname() << "::calibration starting" << std::endl;

  // Get number of layers
  ASSERT(params_.nLayers.value() != boost::none);
  size_t nLayers = *params_.nLayers.value();

  // Allocate data
  weight_.resize(nLayers);
  normalization_.resize(nLayers);
  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    // Create layers
    Layers layers(params_, gdata_, groupName, nx0_, ny0_, groupNz0, nLayers);

    // Save layers
    data_.emplace(groupName, layers);
  }

  // Setup length-scales
  setupLengthScales();

  // Setup weight
  setupWeight();

  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    oops::Log::info() << "Info     : Setup of group " << groupName << ":" << std::endl;

    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      // Setup vertical coordinate
      data_[groupName][jBin].setupVerticalCoord(rv_[groupName], weight_[jBin][groupName]);
    }

    // Print bins
    oops::Log::info() << "Info     : - Horizontal bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[groupName][jBin].rh() << " < ";
    }
    oops::Log::info() << data_[groupName][weight_.size()-1].rh() << std::endl;
    oops::Log::info() << "Info     : - Vertical bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[groupName][jBin].rv() << " < ";
    }
    oops::Log::info() << data_[groupName][weight_.size()-1].rv() << std::endl;

    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      oops::Log::info() << "Info     : * Bin #" << (jBin+1) << ":" << std::endl;

      // Print normalized vertical coordinate
      oops::Log::info() << "Info     :     Normalized vertical coordinate: ";
      for (size_t k0 = 0; k0 < groupNz0-1; ++k0) {
        oops::Log::info() << data_[groupName][jBin].normVertCoord()[k0] << " < ";
      }
      oops::Log::info() << data_[groupName][jBin].normVertCoord()[groupNz0-1]
                        << std::endl;

      // Setup interpolation
      data_[groupName][jBin].setupInterpolation();

      // Setup parallelization
      data_[groupName][jBin].setupParallelization();

      // Setup kernels
      data_[groupName][jBin].setupKernels();

      // Setup normalization
      data_[groupName][jBin].setupNormalization();
      normalization_[jBin].add(data_[groupName][jBin].norm()[groupName]);
    }
  }

  // Get square-root of weight
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    util::sqrtFieldSet(weight_[jBin]);
  }

  oops::Log::trace() << classname() << "::calibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Get number of weight and normalization components from input files
  size_t nWeight = 0;
  size_t nLayers = 0;
  for (size_t ji = 0; ji < inputs_.size(); ++ji) {
    for (size_t jBin = 0; jBin < 100; ++jBin) {
      std::string weightName = "weight_" + std::to_string(jBin);
      if (inputs_[ji].second.name() == weightName) {
        ++nWeight;
      }
      std::string normalizationName = "normalization_" + std::to_string(jBin);
      if (inputs_[ji].second.name() == normalizationName) {
        ++nLayers;
      }
    }
  }

  // Check consistency
  if (nWeight == 0) {
    ASSERT(nLayers == 1);
  } else {
    ASSERT(nLayers == nWeight);
  }

  // Allocate data
  weight_.resize(nLayers);
  normalization_.resize(nLayers);
  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    // Create layers
    Layers layers(params_, gdata_, groupName, nx0_, ny0_, groupNz0, nLayers);

    // Save layers
    data_.emplace(groupName, layers);
  }

  // Get weight and normalization components from input files
  for (size_t ji = 0; ji < inputs_.size(); ++ji) {
    for (size_t jBin = 0; jBin < 100; ++jBin) {
      std::string weightName = "weight_" + std::to_string(jBin);
      if (inputs_[ji].second.name() == weightName) {
        for (const auto & group : groups_) {
          // Group properties
          const std::string groupName = std::get<0>(group);
          const size_t groupNz0 = std::get<1>(group);
          const std::string groupVarInModelFile = std::get<2>(group);

          // Copy field
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(groupName) | atlas::option::levels(groupNz0));
          auto view = atlas::array::make_view<double, 2>(field);
          atlas::Field fieldInput = inputs_[ji].second[groupVarInModelFile];
          auto viewInput = atlas::array::make_view<double, 2>(fieldInput);
          view.assign(viewInput);
          weight_[jBin].add(field);
        }
      }
      std::string normalizationName = "normalization_" + std::to_string(jBin);
      if (inputs_[ji].second.name() == normalizationName) {
        for (const auto & group : groups_) {
          // Group properties
          const std::string groupName = std::get<0>(group);
          const size_t groupNz0 = std::get<1>(group);
          const std::string groupVarInModelFile = std::get<2>(group);

          // Copy field
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(groupName) | atlas::option::levels(groupNz0));
          auto view = atlas::array::make_view<double, 2>(field);
          atlas::Field fieldInput = inputs_[ji].second[groupVarInModelFile];
          auto viewInput = atlas::array::make_view<double, 2>(fieldInput);
          view.assign(viewInput);
          normalization_[jBin].add(field);
        }
      }
    }
  }

  // Set weight to one if not present
  if (nWeight == 0) {
    for (const auto & group : groups_) {
      // Group properties
      const std::string groupName = std::get<0>(group);
      const size_t groupNz0 = std::get<1>(group);

      // Copy field
      atlas::Field field = gdata_.functionSpace().createField<double>(
        atlas::option::name(groupName) | atlas::option::levels(groupNz0));
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(1.0);
      weight_[0].add(field);
    }
  }

  if (comm_.rank() == 0) {
    ASSERT(params_.dataFile.value() != boost::none);

    // NetCDF ids
    int retval, ncid, grpGrpId, layerGrpId;

    // NetCDF file path
    std::string ncfilepath = *params_.dataFile.value();
    ncfilepath.append(".nc");
    oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

    // Open file
    if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

    // Definition mode
    int nLayers;
    if ((retval = nc_get_att_int(ncid, NC_GLOBAL, "nLayers", &nLayers))) ERR(retval);
    ASSERT(nLayers == static_cast<int>(weight_.size()));

    for (const auto & group : groups_) {
      // Group properties
      const std::string groupName = std::get<0>(group);

      // Get group group
      if ((retval = nc_inq_grp_ncid(ncid, groupName.c_str(), &grpGrpId))) ERR(retval);

      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Get layer group
        std::string layerGrpName = "layer_" + std::to_string(jBin);
        if ((retval = nc_inq_grp_ncid(grpGrpId, layerGrpName.c_str(), &layerGrpId))) ERR(retval);
        data_[groupName][jBin].read(layerGrpId);
      }
    }
  }

  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);

    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      // Broadcast data
      data_[groupName][jBin].broadcast();
    }
  }

  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    oops::Log::info() << "Info     : Setup of group " << groupName << ":" << std::endl;

    // Print bins
    oops::Log::info() << "Info     : - Horizontal bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[groupName][jBin].rh() << " < ";
    }
    oops::Log::info() << data_[groupName][weight_.size()-1].rh() << std::endl;
    oops::Log::info() << "Info     : - Vertical bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[groupName][jBin].rv() << " < ";
    }
    oops::Log::info() << data_[groupName][weight_.size()-1].rv() << std::endl;

    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      oops::Log::info() << "Info     : * Bin #" << (jBin+1) << ":" << std::endl;

      // Print normalized vertical coordinate
      oops::Log::info() << "Info     :     Normalized vertical coordinate: ";
      oops::Log::info() << data_[groupName][jBin].normVertCoord()[groupNz0-1]
                        << std::endl;

      // Setup interpolation
      data_[groupName][jBin].setupInterpolation();

      // Setup parallelization
      data_[groupName][jBin].setupParallelization();
    }
  }

  // Get square-root of weight
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    util::sqrtFieldSet(weight_[jBin]);
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (comm_.rank() == 0 && (params_.dataFile.value() != boost::none)) {
    // NetCDF ids
    int retval, ncid, grpGrpId, layerGrpId;
    std::vector<std::vector<int>> grpIdsVec;

    // NetCDF file path
    std::string ncfilepath = *params_.dataFile.value();
    ncfilepath.append(".nc");
    oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

    // Create file
    if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid))) ERR(retval);

    // Definition mode
    int nLayers = weight_.size();
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "nLayers", NC_INT, 1, &nLayers))) ERR(retval);

    for (const auto & group : groups_) {
      // Group properties
      const std::string groupName = std::get<0>(group);

      // Create group group
      if ((retval = nc_def_grp(ncid, groupName.c_str(), &grpGrpId))) ERR(retval);

      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create layer group
        std::string layerGrpName = "layer_" + std::to_string(jBin);
        if ((retval = nc_def_grp(grpGrpId, layerGrpName.c_str(), &layerGrpId))) ERR(retval);
        std::vector<int> grpIds = data_.at(groupName)[jBin].writeDef(layerGrpId);
        grpIdsVec.push_back(grpIds);
      }
    }

    // End definition mode
    if ((retval = nc_enddef(ncid))) ERR(retval);

    // Data mode
    size_t jj = 0;
    for (const auto & group : groups_) {
      // Group properties
      const std::string groupName = std::get<0>(group);

      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create layer group
        data_.at(groupName)[jBin].writeData(grpIdsVec[jj]);
        ++jj;
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> FastLAM::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Create vector of pairs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> pairs;

  // Get vector of configurations
  std::vector<eckit::LocalConfiguration> outputModelFilesConf
    = params_.outputModelFilesConf.value().get_value_or({});

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  for (const auto & conf : outputModelFilesConf) {
    // Get file configuration
    eckit::LocalConfiguration file = getFileConf(comm_, conf);

    // Get parameter
    const std::string param = conf.getString("parameter");

    if (param == "index i") {
      // Copy field
      atlas::FieldSet fset;
      atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
      auto viewIndexI0 = atlas::array::make_view<int, 1>(fieldIndexI0);
      for (const auto & var : activeVars_.variables()) {
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(var) | atlas::option::levels(activeVars_.getLevels(var)));
        auto view = atlas::array::make_view<double, 2>(field);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (viewGhost0(jnode0) == 0) {
            for (int k0 = 0; k0 < activeVars_.getLevels(var); ++k0) {
              view(jnode0, k0) = static_cast<double>(viewIndexI0(jnode0));
            }
          }
        }
        fset.add(field);
      }

      // Add pair
      pairs.push_back(std::make_pair(file, fset));
    }
    if (param == "index j") {
      // Copy field
      atlas::FieldSet fset;
      atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
      auto viewIndexJ0 = atlas::array::make_view<int, 1>(fieldIndexJ0);
      for (const auto & var : activeVars_.variables()) {
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(var) | atlas::option::levels(activeVars_.getLevels(var)));
        auto view = atlas::array::make_view<double, 2>(field);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (viewGhost0(jnode0) == 0) {
            for (int k0 = 0; k0 < activeVars_.getLevels(var); ++k0) {
              view(jnode0, k0) = static_cast<double>(viewIndexJ0(jnode0));
            }
          }
        }
        fset.add(field);
      }

      // Add pair
      pairs.push_back(std::make_pair(file, fset));
    }
    if (param == "normalized horizontal length-scale") {
      // Copy fields
      atlas::FieldSet fset;
      for (const auto & var : activeVars_.variables()) {
        // Default: missing value
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(var) | atlas::option::levels(activeVars_.getLevels(var)));
        auto view = atlas::array::make_view<double, 2>(field);
        view.assign(util::missingValue<double>());
        fset.add(field);

        for (const auto & group : groups_) {
          // Group properties
          const std::string groupName = std::get<0>(group);
          const std::string groupVarInModelFile = std::get<2>(group);

          // Copy field
          if (groupVarInModelFile == var) {
            const atlas::Field fieldRh = rh_[groupName];
            const auto viewRh = atlas::array::make_view<double, 2>(fieldRh);
            view.assign(viewRh);
          }
        }
      }

      // Add pair
      pairs.push_back(std::make_pair(file, fset));
    }
    if (param == "weight") {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Copy fields
        atlas::FieldSet fset;
        for (const auto & var : activeVars_.variables()) {
          // Default: missing value
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(var) | atlas::option::levels(activeVars_.getLevels(var)));
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(util::missingValue<double>());
          fset.add(field);

          for (const auto & group : groups_) {
            // Group properties
            const std::string groupName = std::get<0>(group);
            const std::string groupVarInModelFile = std::get<2>(group);

            // Copy field
            if (groupVarInModelFile == var) {
              const atlas::Field fieldWgtSqrt = weight_[jBin][groupName];
              const auto viewWgtSqrt = atlas::array::make_view<double, 2>(fieldWgtSqrt);
              for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
                if (viewGhost0(jnode0) == 0) {
                  for (int k0 = 0; k0 < activeVars_.getLevels(var); ++k0) {
                    view(jnode0, k0) = viewWgtSqrt(jnode0, k0)*viewWgtSqrt(jnode0, k0);
                  }
                }
              }
            }
          }
        }

        // Update configuration
        eckit::LocalConfiguration fileCmp(file);
        util::seekAndReplace(fileCmp, "%component%", jBin, 2);

        // Add pair
        pairs.push_back(std::make_pair(fileCmp, fset));
      }
    }
    if (param == "normalization") {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Copy fields
        atlas::FieldSet fset;
        for (const auto & var : activeVars_.variables()) {
          // Default: missing value
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(var) | atlas::option::levels(activeVars_.getLevels(var)));
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(util::missingValue<double>());
          fset.add(field);

          for (const auto & group : groups_) {
            // Group properties
            const std::string groupName = std::get<0>(group);
            const std::string groupVarInModelFile = std::get<2>(group);

            // Copy field
            if (groupVarInModelFile == var) {
              const atlas::Field fieldNorm = normalization_[jBin][groupName];
              const auto viewNorm = atlas::array::make_view<double, 2>(fieldNorm);
              view.assign(viewNorm);
            }
          }
        }

        // Update configuration
        eckit::LocalConfiguration fileCmp(file);
        util::seekAndReplace(fileCmp, "%component%", jBin, 2);

        // Add pair
        pairs.push_back(std::make_pair(fileCmp, fset));
      }
    }
    if (param == "normalization accuracy") {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Copy fields
        atlas::FieldSet fset;
        for (const auto & var : activeVars_.variables()) {
          // Default: missing value
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(var) | atlas::option::levels(activeVars_.getLevels(var)));
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(util::missingValue<double>());
          fset.add(field);

          for (const auto & group : groups_) {
            // Group properties
            const std::string groupName = std::get<0>(group);
            const std::string groupVarInModelFile = std::get<2>(group);

            // Copy field
            if (groupVarInModelFile == var) {
              const atlas::Field fieldNorm = data_.at(groupName)[jBin].normAcc()[groupName];
              const auto viewNorm = atlas::array::make_view<double, 2>(fieldNorm);
              view.assign(viewNorm);
            }
          }
        }

        // Update configuration
        eckit::LocalConfiguration fileCmp(file);
        util::seekAndReplace(fileCmp, "%component%", jBin, 2);

        // Add pair
        pairs.push_back(std::make_pair(fileCmp, fset));
      }
    }
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return pairs;
}

// -----------------------------------------------------------------------------

void FastLAM::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

void FastLAM::setupLengthScales() {
  oops::Log::trace() << classname() << "::setupLengthScales starting" << std::endl;

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Get rh and rv from input files
  for (size_t ji = 0; ji < inputs_.size(); ++ji) {
    if (inputs_[ji].second.name() == "rh") {
      for (const auto & group : groups_) {
        // Group properties
        const std::string groupName = std::get<0>(group);
        const size_t groupNz0 = std::get<1>(group);
        const std::string groupVarInModelFile = std::get<2>(group);

        // Copy field
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(groupName) | atlas::option::levels(groupNz0));
        auto view = atlas::array::make_view<double, 2>(field);
        if (inputs_[ji].second.has_field(groupVarInModelFile)) {
          const atlas::Field fieldInput = inputs_[ji].second[groupVarInModelFile];
          const auto viewInput = atlas::array::make_view<double, 2>(fieldInput);
          view.assign(viewInput);
        } else {
          throw eckit::Exception("missing " + groupVarInModelFile + " in model file", Here());
        }
        rh_.add(field);
      }
    }
    if (inputs_[ji].second.name() == "rv") {
      for (const auto & group : groups_) {
        // Group properties
        const std::string groupName = std::get<0>(group);
        const size_t groupNz0 = std::get<1>(group);
        const std::string groupVarInModelFile = std::get<2>(group);

        // Copy field
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(groupName) | atlas::option::levels(groupNz0));
        auto view = atlas::array::make_view<double, 2>(field);
        if (inputs_[ji].second.has_field(groupVarInModelFile)) {
          const atlas::Field fieldInput = inputs_[ji].second[groupVarInModelFile];
          const auto viewInput = atlas::array::make_view<double, 2>(fieldInput);
          view.assign(viewInput);
        } else {
          throw eckit::Exception("missing " + groupVarInModelFile + " in model file", Here());
        }
        rv_.add(field);
      }
    }
  }

  // Get rh and rv from yaml
  if (rh_.empty()) {
    ASSERT(params_.rhFromYaml.value() != boost::none);
    for (const auto & group : groups_) {
      // Group properties
      const std::string groupName = std::get<0>(group);
      const size_t groupNz0 = std::get<1>(group);

      // Get yaml value/profile
      std::vector<double> profile;
      for (const auto & vParams : *params_.rhFromYaml.value()) {
        const std::string vGroupName = vParams.group.value();
        if (vGroupName == groupName) {
          if (profile.size() > 0) {
            throw eckit::UserError("group" + groupName + " present more that once", Here());
          }
          profile.resize(groupNz0);
          if (vParams.value.value() != boost::none
            && vParams.profile.value() != boost::none) {
            throw eckit::UserError("both value and profile present in the same item", Here());
          }
          if (vParams.value.value() != boost::none) {
            // Copy value
            std::fill(profile.begin(), profile.end(), *vParams.value.value());
          } else if (vParams.profile.value() != boost::none) {
            // Copy profile
            ASSERT(vParams.profile.value()->size() == groupNz0);
            profile = *vParams.profile.value();
          } else {
            throw eckit::UserError("missing value or profile for group " + groupName, Here());
          }
        }
      }
      ASSERT(profile.size() > 0);

      // Copy to rh_
      atlas::Field fieldRh = gdata_.functionSpace().createField<double>(
        atlas::option::name(groupName) | atlas::option::levels(groupNz0));
      auto viewRh = atlas::array::make_view<double, 2>(fieldRh);
      for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
        if (viewGhost0(jnode0) == 0) {
          for (size_t k0 = 0; k0 < groupNz0; ++k0) {
            viewRh(jnode0, k0) = profile[k0];
          }
        }
      }
      rh_.add(fieldRh);
    }
  }

  if (rv_.empty()) {
    ASSERT(params_.rvFromYaml.value() != boost::none);
    for (const auto & group : groups_) {
      // Group properties
      const std::string groupName = std::get<0>(group);
      const size_t groupNz0 = std::get<1>(group);

      // Get yaml value/profile
      std::vector<double> profile;
      for (const auto & vParams : *params_.rvFromYaml.value()) {
        const std::string vGroupName = vParams.group.value();
        if (vGroupName == groupName) {
          if (profile.size() > 0) {
            throw eckit::UserError("group" + groupName + " present more that once", Here());
          }
          profile.resize(groupNz0);
          if (vParams.value.value() != boost::none
            && vParams.profile.value() != boost::none) {
            throw eckit::UserError("both value and profile present in the same item", Here());
          }
          if (vParams.value.value() != boost::none) {
            // Copy value
            std::fill(profile.begin(), profile.end(), *vParams.value.value());
          } else if (vParams.profile.value() != boost::none) {
            // Copy profile
            ASSERT(vParams.profile.value()->size() == groupNz0);
            profile = *vParams.profile.value();
          } else {
            throw eckit::UserError("missing value or profile for group " + groupName, Here());
          }
        }
      }
      ASSERT(profile.size() > 0);

      // Copy to rv_
      atlas::Field fieldRv = gdata_.functionSpace().createField<double>(
        atlas::option::name(groupName) | atlas::option::levels(groupNz0));
      auto viewRv = atlas::array::make_view<double, 2>(fieldRv);
      for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
        if (viewGhost0(jnode0) == 0) {
          for (size_t k0 = 0; k0 < groupNz0; ++k0) {
            viewRv(jnode0, k0) = profile[k0];
          }
        }
      }
      rv_.add(fieldRv);
    }
  }

  // Compute correlation
  oops::Log::info() << "Info     : Correlation rh / rv:" << std::endl;
  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    // Get fields
    atlas::Field fieldRh = rh_[groupName];
    atlas::Field fieldRv = rv_[groupName];
    auto viewRh = atlas::array::make_view<double, 2>(fieldRh);
    auto viewRv = atlas::array::make_view<double, 2>(fieldRv);

    // Compute mean
    double meanRh = 0.0;
    double meanRv = 0.0;
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        for (size_t k0 = 0; k0 < groupNz0; ++k0) {
          meanRh += viewRh(jnode0, k0);
          meanRv += viewRv(jnode0, k0);
        }
      }
    }
    comm_.allReduceInPlace(meanRh, eckit::mpi::sum());
    comm_.allReduceInPlace(meanRv, eckit::mpi::sum());
    meanRh *= 1.0/static_cast<double>(nx0_*ny0_*groupNz0);
    meanRv *= 1.0/static_cast<double>(nx0_*ny0_*groupNz0);

    // Compute horizontal variance and covariance
    double varRh = 0.0;
    double varRv = 0.0;
    double cov = 0.0;
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        for (size_t k0 = 0; k0 < groupNz0; ++k0) {
          varRh += (viewRh(jnode0, k0)-meanRh)*(viewRh(jnode0, k0)-meanRh);
          varRv += (viewRv(jnode0, k0)-meanRv)*(viewRv(jnode0, k0)-meanRv);
          cov += (viewRh(jnode0, k0)-meanRh)*(viewRv(jnode0, k0)-meanRv);
        }
      }
    }
    comm_.allReduceInPlace(varRh, eckit::mpi::sum());
    comm_.allReduceInPlace(varRv, eckit::mpi::sum());
    comm_.allReduceInPlace(cov, eckit::mpi::sum());

    // Compute horizontal correlation
    double cor = 0.0;
    if (varRh > 0.0 && varRv > 0.0) {
      cor = cov/std::sqrt(varRh*varRv);
    }

    // Print result
    oops::Log::info() << "Info     : - Group " << groupName << ": " << cor << std::endl;
  }

  oops::Log::trace() << classname() << "::setupLengthScales done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::setupWeight() {
  oops::Log::trace() << classname() << "::setupWeight starting" << std::endl;

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    // Normalize rh
    atlas::Field fieldRh = rh_[groupName];
    auto viewRh = atlas::array::make_view<double, 2>(fieldRh);
    auto viewArea =
      atlas::array::make_view<double, 2>(gdata_.fieldSet()["area"]);
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        for (size_t k0 = 0; k0 < groupNz0; ++k0) {
          viewRh(jnode0, k0) = viewRh(jnode0, k0)/std::sqrt(viewArea(jnode0, 0));
        }
      }
    }
  }

  // Initialize weight
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    atlas::FieldSet fsetWeight = util::copyFieldSet(rh_);
    util::zeroFieldSet(fsetWeight);
    weight_[jBin] = fsetWeight;
  }

  for (const auto & group : groups_) {
    // Group properties
    const std::string groupName = std::get<0>(group);
    const size_t groupNz0 = std::get<1>(group);

    // Get rh field general properties
    atlas::Field fieldRh = rh_[groupName];
    auto viewRh = atlas::array::make_view<double, 2>(fieldRh);
    double minRh = std::numeric_limits<double>::max();
    double maxRh = std::numeric_limits<double>::min();
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        for (size_t k0 = 0; k0 < groupNz0; ++k0) {
          minRh = std::min(minRh, viewRh(jnode0, k0));
          maxRh = std::max(maxRh, viewRh(jnode0, k0));
        }
      }
    }
    comm_.allReduceInPlace(minRh, eckit::mpi::min());
    comm_.allReduceInPlace(maxRh, eckit::mpi::max());

    if (weight_.size() > 1) {
      // Prepare binning
      double binWidth = (maxRh-minRh)/static_cast<double>(weight_.size()-1);
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        data_[groupName][jBin].rh() = minRh+static_cast<double>(jBin)*binWidth;
      }

      // Compute weight
      for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
        if (viewGhost0(jnode0) == 0) {
          for (size_t k0 = 0; k0 < groupNz0; ++k0) {
            // Raw weight (difference-based)
            std::vector<double> weight(weight_.size());
            for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
              atlas::Field fieldWgt = weight_[jBin][groupName];
              auto viewWgt = atlas::array::make_view<double, 2>(fieldWgt);
              double diff = std::abs(viewRh(jnode0, k0)-data_[groupName][jBin].rh())/(maxRh-minRh);
              viewWgt(jnode0, k0) = std::exp(-4.6*diff);  // Factor 4.6 => minimum weight ~0.01
            }

            // Normalized weight
            double weightSum = 0.0;
            for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
              atlas::Field fieldWgt = weight_[jBin][groupName];
              const auto viewWgt = atlas::array::make_view<double, 2>(fieldWgt);
              weightSum += viewWgt(jnode0, k0);
            }
            for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
              atlas::Field fieldWgt = weight_[jBin][groupName];
              auto viewWgt = atlas::array::make_view<double, 2>(fieldWgt);
              viewWgt(jnode0, k0) = viewWgt(jnode0, k0)/weightSum;
            }
          }
        }
      }
    } else {
      // Single bin
      data_[groupName][0].rh() = 0.5*(minRh+maxRh);

      // Set weight to one
      atlas::Field fieldWgt = weight_[0][groupName];
      auto viewWgt = atlas::array::make_view<double, 2>(fieldWgt);
      viewWgt.assign(1.0);
    }
  }

  oops::Log::trace() << classname() << "::setupWeight done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t FastLAM::getGroupIndex(const std::string & var) const {
  oops::Log::trace() << classname() << "::getGroupIndex starting" << std::endl;

  // Initialization
  size_t groupIndex = 0;

  for (const auto & group : groups_) {
    // Group properties
    const std::vector<std::string> groupVars = std::get<3>(group);

    // Check group variables
    if (std::find(groupVars.begin(), groupVars.end(), var) != groupVars.end()) {
      break;
    }

    // Update group index
    ++groupIndex;
  }

  // Check group index has been found
  ASSERT(groupIndex < groups_.size());

  oops::Log::trace() << classname() << "::getGroupIndex done" << std::endl;
  return groupIndex;
}

// -----------------------------------------------------------------------------

size_t FastLAM::getK0Offset(const std::string & var) const {
  oops::Log::trace() << classname() << "::getK0Offset starting" << std::endl;

  // Default value
  size_t k0Offset = 0;

  if (std::find(active2dVars_.variables().begin(), active2dVars_.variables().end(), var)
    != active2dVars_.variables().end()) {
    // This is a 2d variable
    if (params_.lev2d.value() == "last") {
      // Use the last level of the group
      k0Offset = std::get<1>(groups_[getGroupIndex(var)])-1;
    }
  }

  oops::Log::trace() << classname() << "::getK0Offset starting" << std::endl;
  return k0Offset;
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration FastLAM::getFileConf(const eckit::mpi::Comm & comm,
                                               const eckit::Configuration & conf) const {
  oops::Log::trace() << classname() << "::getFileConf starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm.size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif

  // Get IO configuration
  eckit::LocalConfiguration file = conf.getSubConfiguration("file");

  // Replace patterns
  util::seekAndReplace(file, "_MPI_", mpi);
  util::seekAndReplace(file, "_OMP_", omp);

  oops::Log::trace() << classname() << "::getFileConf done" << std::endl;
  return file;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber



