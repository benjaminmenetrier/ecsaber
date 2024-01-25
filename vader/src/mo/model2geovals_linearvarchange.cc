/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"

#include "mo/constants.h"
#include "mo/model2geovals_linearvarchange.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;


namespace mo {

void evalAirPressureTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds) {
  // height at pressure levels
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  // height at theta levels
  const auto hView = make_view<const double, 2>(augStateFlds["height"]);

  const auto pIncView = make_view<double, 2>(incFlds["air_pressure_levels"]);
  auto pbarIncView = make_view<double, 2>(incFlds["air_pressure"]);

  double alpha_jl(0.0);

  atlas::idx_t lvls(incFlds["air_pressure"].levels());
  for (atlas::idx_t jn = 0; jn < pbarIncView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      alpha_jl = (hView(jn, jl) - hlView(jn, jl)) / (hlView(jn, jl+1) - hlView(jn, jl));
      pbarIncView(jn, jl) = (1.0 - alpha_jl) * pIncView(jn, jl) +
                            alpha_jl * pIncView(jn, jl+1);
    }
  }
}

void evalAirPressureAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds) {
  // height at pressure levels
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  // height at theta levels
  const auto hView = make_view<const double, 2>(augStateFlds["height"]);

  auto pbarHatView = make_view<double, 2>(hatFlds["air_pressure"]);
  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);

  double alpha_jl(0.0);

  atlas::idx_t lvls(hatFlds["air_pressure"].levels());
  for (atlas::idx_t jn = 0; jn < pHatView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      alpha_jl = (hView(jn, jl) - hlView(jn, jl)) / (hlView(jn, jl+1) - hlView(jn, jl));
      pHatView(jn, jl) += (1.0 - alpha_jl) * pbarHatView(jn, jl);
      pHatView(jn, jl+1) += alpha_jl * pbarHatView(jn, jl);
      pbarHatView(jn, jl) = 0.0;
    }
  }
}

}  // namespace mo
