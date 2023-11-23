/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/constants.h"
#include "mo/eval_exner.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {

// ------------------------------------------------------------------------------------------------
// Note this is a copy of evalHydrostaticExner
void eval_hydrostatic_exner_levels_tl(atlas::FieldSet & incFlds,
                                      const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_hydrostatic_exner_levels_tl()] starting ..." << std::endl;
  const auto pView = make_view<const double, 2>(augStateFlds["hydrostatic_pressure_levels"]);
  const auto exnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  const auto pIncView = make_view<const double, 2>(incFlds["hydrostatic_pressure_levels"]);
  auto exnerIncView = make_view<double, 2>(incFlds["hydrostatic_exner_levels"]);

  atlas::idx_t levels = incFlds["hydrostatic_exner_levels"].levels();
  for (atlas::idx_t jn = 0; jn < incFlds["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < levels; ++jl) {
      exnerIncView(jn, jl) = pIncView(jn, jl) *
        (constants::rd_over_cp * exnerView(jn, jl)) /
        pView(jn, jl);
    }
  }
  oops::Log::trace() << "[eval_hydrostatic_exner_levels_tl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_hydrostatic_exner_levels_ad(atlas::FieldSet & hatFlds,
                                      const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_hydrostatic_exner_levels_ad()] starting ..." << std::endl;
  const auto pView = make_view<const double, 2>(augStateFlds["hydrostatic_pressure_levels"]);
  const auto exnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  auto pHatView = make_view<double, 2>(hatFlds["hydrostatic_pressure_levels"]);
  auto exnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  atlas::idx_t levels = hatFlds["hydrostatic_exner_levels"].levels();
  for (atlas::idx_t jn = 0; jn < hatFlds["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < levels; ++jl) {
      pHatView(jn, jl) += exnerHatView(jn, jl) *
        (constants::rd_over_cp * exnerView(jn, jl)) /
        pView(jn, jl);
      exnerHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_hydrostatic_exner_levels_ad()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_hydrostatic_exner_levels_tl_inv(atlas::FieldSet & incFlds,
                                          const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_hydrostatic_exner_levels_tl_inv()] starting ..." << std::endl;
  const auto pView = make_view<const double, 2>(augStateFlds["hydrostatic_pressure_levels"]);
  const auto exnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  const auto exnerIncView = make_view<const double, 2>(incFlds["hydrostatic_exner_levels"]);
  auto pIncView = make_view<double, 2>(incFlds["hydrostatic_pressure_levels"]);

  for (atlas::idx_t jn = 0; jn < pIncView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < pIncView.shape(1); ++jl) {
      pIncView(jn, jl) = exnerIncView(jn, jl) * pView(jn, jl) /
              (constants::rd_over_cp * exnerView(jn, jl));
    }
  }
  oops::Log::trace() << "[eval_hydrostatic_exner_levels_tl_inv()] ... done" << std::endl;
}

}  // namespace mo
