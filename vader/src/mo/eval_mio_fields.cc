/*
 * (C) Crown Copyright 2023 Met Office
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
#include "mo/eval_mio_fields.h"
#include "mo/functions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {

void eval_mio_fields_nl(const std::string & mio_path, atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_mio_fields_nl()] starting ..." << std::endl;
  const auto rhtView = make_view<const double, 2>(augStateFlds["rht"]);
  const auto clView = make_view<const double, 2>
                (augStateFlds["liquid_cloud_volume_fraction_in_atmosphere_layer"]);
  const auto cfView = make_view<const double, 2>
                (augStateFlds["ice_cloud_volume_fraction_in_atmosphere_layer"]);

  auto cleffView = make_view<double, 2>(augStateFlds["cleff"]);
  auto cfeffView = make_view<double, 2>(augStateFlds["cfeff"]);

  Eigen::MatrixXd mioCoeffCl = functions::createMIOCoeff(mio_path, "qcl_coef");
  Eigen::MatrixXd mioCoeffCf = functions::createMIOCoeff(mio_path, "qcf_coef");

  for  (atlas::idx_t jn = 0; jn < augStateFlds["rht"].shape(0); ++jn) {
    for (int jl = 0; jl < augStateFlds["rht"].levels(); ++jl) {
      if (jl < constants::mioLevs) {
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
  oops::Log::trace() << "[eval_mio_fields_nl()] ... exit" << std::endl;
}

}  // namespace mo
