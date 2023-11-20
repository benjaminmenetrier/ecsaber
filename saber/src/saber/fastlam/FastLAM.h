/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/fastlam/FastLAMParametersBase.h"
#include "saber/fastlam/Layers.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

class FastLAMParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(FastLAMParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<FastLAMParametersBase> read{"read", this};
  oops::OptionalParameter<FastLAMParametersBase> calibration{"calibration", this};

  oops::patch::Variables mandatoryActiveVars() const override {return oops::patch::Variables();}
};

// -----------------------------------------------------------------------------

class FastLAM : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::fastlam::FastLAM";}

  typedef FastLAMParameters     Parameters_;
  typedef FastLAMParametersBase ParametersBase_;

  FastLAM(const oops::GeometryData &,
          const oops::patch::Variables &,
          const eckit::Configuration &,
          const Parameters_ &,
          const oops::FieldSet3D &,
          const oops::FieldSet3D &);

  virtual ~FastLAM();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToRead() override;

  void read() override;

  void directCalibration(const std::vector<atlas::FieldSet> &) override;

  void write() const override;
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite() const override;

 private:
  // Model grid geometry data
  const oops::GeometryData & gdata_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Active variables
  const oops::patch::Variables activeVars_;

  // Active 2D variables
  oops::patch::Variables active2dVars_;

  // Groups of variables
  std::vector<std::tuple<std::string, size_t, std::string, std::vector<std::string>>> groups_;

  // Parameters
  ParametersBase_ params_;

  // Inputs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> inputs_;
  atlas::FieldSet rh_;
  atlas::FieldSet rv_;

  // Weight and normalization
  std::vector<atlas::FieldSet> weight_;
  std::vector<atlas::FieldSet> normalization_;

  // Data
  std::unordered_map<std::string, Layers> data_;

  // Model grid
  size_t nx0_;
  size_t ny0_;
  size_t nodes0_;

  // Setup length-scales
  void setupLengthScales();

  // Setup weight
  void setupWeight();

  // Utilities
  size_t getGroupIndex(const std::string &) const;
  size_t getK0Offset(const std::string &) const;
  eckit::LocalConfiguration getFileConf(const eckit::mpi::Comm &,
                                        const eckit::Configuration &) const;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber