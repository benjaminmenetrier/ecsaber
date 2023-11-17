/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "mo/constants.h"
#include "mo/eval_air_pressure_levels.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

/// \details Vertical extrapolation of air pressure above model top
void eval_air_pressure_levels_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_air_pressure_levels_nl()] starting ..." << std::endl;

  auto ds_elmo = make_view<const double, 2>(stateFlds["exner_levels_minus_one"]);
  auto ds_plmo = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto ds_t = make_view<const double, 2>(stateFlds["potential_temperature"]);
  auto ds_hl = make_view<const double, 2>(stateFlds["height_levels"]);
  auto ds_pl = make_view<double, 2>(stateFlds["air_pressure_levels"]);

  idx_t levels(stateFlds["air_pressure_levels"].levels());
  for (idx_t jn = 0; jn < stateFlds["air_pressure_levels"].shape(0); ++jn) {
    for (idx_t jl = 0; jl < levels - 1; ++jl) {
      ds_pl(jn, jl) = ds_plmo(jn, jl);
    }

    // Note that I am calculating the exner pressure above the top first and then
    // converting it to pressure
    // Note that strictly speaking we should be using virtual potential temperature here
    // but given that there should be no moisture at the top of the model, we should be
    // able to ignore that contribution.
    //
    // exner^k+1 = exner^k - g (height^k+1 - height^k) / (cp theta_v)
    //
    // pressure^k+1 = reference_pressure * (exner^k+1)**((1.0 / constants::rd_over_cp)
    //
    // where k is the model level index on half levels just below model top.

    ds_pl(jn, levels-1) =  constants::p_zero * pow(
      ds_elmo(jn, levels-2) - (constants::grav * (ds_hl(jn, levels-1) - ds_hl(jn, levels-2))) /
      (constants::cp * ds_t(jn, levels-2)), (1.0 / constants::rd_over_cp));

    ds_pl(jn, levels-1) = ds_pl(jn, levels-1) > 0.0 ? ds_pl(jn, levels-1) : constants::deps;
  }

  oops::Log::trace() << "[eval_air_pressure_levels_nl()] ... exit" << std::endl;
}

void eval_air_pressure_levels_tl(atlas::FieldSet & incFlds, const atlas::FieldSet & stateFlds) {
  eckit::NotImplemented(Here());
}

void eval_air_pressure_levels_ad(atlas::FieldSet & incFlds, const atlas::FieldSet & stateFlds) {
  eckit::NotImplemented(Here());
}

}  // namespace mo
