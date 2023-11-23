/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Core>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/constants.h"
#include "mo/functions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {
namespace functions {

std::vector<double> getLookUp(const std::string & sVPFilePath,
                              const std::string & shortName,
                              const std::size_t lookupSize) {
  std::vector<double> values(lookupSize, 0);

  umGetLookUp_f90(static_cast<int>(sVPFilePath.size()),
                  sVPFilePath.c_str(),
                  static_cast<int>(shortName.size()),
                  shortName.c_str(),
                  static_cast<int>(lookupSize),
                  values[0]);

  return std::vector<double>(values.begin(), values.end());
}


std::vector<std::vector<double>> getLookUps(const std::string & sVPFilePath,
                                            const oops::patch::Variables & vars,
                                            const std::size_t lookupSize) {
  std::vector<std::vector<double>> values(vars.size(), std::vector<double>(lookupSize));

  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = getLookUp(sVPFilePath, vars[i], lookupSize);
  }

  return values;
}

void getMIOFields(atlas::FieldSet & augStateFlds) {
  const auto rhtView = make_view<const double, 2>(augStateFlds["rht"]);
  const auto clView = make_view<const double, 2>
                (augStateFlds["liquid_cloud_volume_fraction_in_atmosphere_layer"]);
  const auto cfView = make_view<const double, 2>
                (augStateFlds["ice_cloud_volume_fraction_in_atmosphere_layer"]);

  auto cleffView = make_view<double, 2>(augStateFlds["cleff"]);
  auto cfeffView = make_view<double, 2>(augStateFlds["cfeff"]);

  Eigen::MatrixXd mioCoeffCl = createMIOCoeff(constants::mioCoefficientsFilePath, "qcl_coef");
  Eigen::MatrixXd mioCoeffCf = createMIOCoeff(constants::mioCoefficientsFilePath, "qcf_coef");

  for  (atlas::idx_t jn = 0; jn < augStateFlds["rht"].shape(0); ++jn) {
    for (int jl = 0; jl < augStateFlds["rht"].levels(); ++jl) {
      if (jl < static_cast<int>(constants::mioLevs)) {
        std::size_t ibin = (rhtView(jn, jl) > 1.0) ? constants::mioBins - 1 :
                           static_cast<std::size_t>(floor(rhtView(jn, jl) / constants::rHTBin));

        std::double_t ceffdenom = (1.0 -  clView(jn, jl) * cfView(jn, jl) );
        if (ceffdenom > constants::tol) {
          std::double_t clcf = clView(jn, jl) * cfView(jn, jl);
          cleffView(jn, jl) = mioCoeffCl(jl, ibin) * (clView(jn, jl) - clcf) / ceffdenom;
          cfeffView(jn, jl) = mioCoeffCf(jl, ibin) * (cfView(jn, jl) - clcf) / ceffdenom;
        } else {
          cleffView(jn, jl) = 0.5;
          cfeffView(jn, jl) = 0.5;
        }
      } else {
        cleffView(jn, jl) = 0.0;
        cfeffView(jn, jl) = 0.0;
      }
    }
  }
}

Eigen::MatrixXd createMIOCoeff(const std::string mioFileName,
                               const std::string s) {
  oops::Log::trace() << "[createMIOCoeff] starting ..." << std::endl;
  Eigen::MatrixXd mioCoeff(static_cast<std::size_t>(constants::mioLevs),
                           static_cast<std::size_t>(constants::mioBins));

  std::vector<double> valuesvec(constants::mioLookUpLength, 0);

  umGetLookUp2D_f90(static_cast<int>(mioFileName.size()),
                    mioFileName.c_str(),
                    static_cast<int>(s.size()),
                    s.c_str(),
                    static_cast<int>(constants::mioBins),
                    static_cast<int>(constants::mioLevs),
                    valuesvec[0]);

  for (int j = 0; j < static_cast<int>(constants::mioLevs); ++j) {
    for (int i = 0; i < static_cast<int>(constants::mioBins); ++i) {
      // Fortran returns column major order, but C++ needs row major
      mioCoeff(j, i) = valuesvec[i * constants::mioLevs+j];
    }
  }
  oops::Log::trace() << "[createMIOCoeff] ... exit" << std::endl;
  return mioCoeff;
}

}  // namespace functions
}  // namespace mo
