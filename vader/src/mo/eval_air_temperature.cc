/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "mo/constants.h"
#include "mo/eval_air_temperature.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

/// \details Calculates air temperature from potential temperature and exner.
void eval_air_temperature_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_air_temperature_nl()] starting ..." << std::endl;

  const auto ds_theta = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto ds_exner = make_view<const double, 2>(stateFlds["exner"]);
  auto ds_atemp = make_view<double, 2>(stateFlds["air_temperature"]);

  auto fspace = stateFlds["air_temperature"].functionspace();
  idx_t lvls(stateFlds["air_temperature"].levels());
  for (atlas::idx_t jn = 0; jn < ds_atemp.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      ds_atemp(jn, jl) = ds_theta(jn, jl) * ds_exner(jn, jl);
    }
  }

  oops::Log::trace() << "[eval_air_temperature_nl()] ... exit" << std::endl;
}

void eval_air_temperature_tl(atlas::FieldSet & incFlds, const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_air_temperature_tl()] starting ..." << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto exnerLevelsView = make_view<const double, 2>(stateFlds["exner_levels_minus_one"]);
  const auto thetaView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  const auto exnerLevelsIncView = make_view<const double, 2>(incFlds["exner_levels_minus_one"]);
  const auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  auto tIncView = make_view<double, 2>(incFlds["air_temperature"]);

  idx_t lvls(incFlds["air_temperature"].levels());
  idx_t lvlsm1 = lvls - 1;
  double exnerTopVal;
  double exnerTopIncVal;

  // Active code
  for (idx_t jn = 0; jn < tIncView.shape(0); ++jn) {
    // Passive code: Value above model top is assumed to be in hydrostatic balance.
    exnerTopVal = exnerLevelsView(jn, lvlsm1) -
      (constants::grav * (hlView(jn, lvls) - hlView(jn, lvlsm1))) /
      (constants::cp * thetaView(jn, lvlsm1));

    // Active code;
    for (idx_t jl= 0; jl < lvls - 1; ++jl) {
       tIncView(jn, jl) = (
         ( (hView(jn, jl) - hlView(jn, jl)) * exnerLevelsView(jn, jl + 1) +
           (hlView(jn, jl+1)  - hView(jn, jl)) * exnerLevelsView(jn, jl) ) *
           thetaIncView(jn, jl) +
         ( (hView(jn, jl) - hlView(jn, jl)) * exnerLevelsIncView(jn, jl + 1) +
           (hlView(jn, jl+1)  - hView(jn, jl)) * exnerLevelsIncView(jn, jl) ) *
         thetaView(jn, jl) ) /
         (hlView(jn, jl+1) - hlView(jn, jl));
    }

    exnerTopIncVal = exnerLevelsIncView(jn, lvlsm1) +
      thetaIncView(jn, lvlsm1) * (exnerLevelsView(jn, lvlsm1) - exnerTopVal) /
      thetaView(jn, lvlsm1);

    tIncView(jn, lvlsm1) = (
      ( (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) * exnerTopVal +
        (hlView(jn, lvls)  - hView(jn, lvlsm1)) * exnerLevelsView(jn, lvlsm1) ) *
         thetaIncView(jn, lvlsm1) +
      ( (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) * exnerTopIncVal +
        (hlView(jn, lvls)  - hView(jn, lvlsm1)) * exnerLevelsIncView(jn, lvlsm1) ) *
         thetaView(jn, lvlsm1)) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));
  }
  oops::Log::trace() << "[eval_air_temperature_tl()] ... exit" << std::endl;
}

void eval_air_temperature_ad(atlas::FieldSet & hatFlds, const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_air_temperature_ad()] starting ..." << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto exnerLevelsView = make_view<const double, 2>(stateFlds["exner_levels_minus_one"]);
  const auto thetaView = make_view<const double, 2>(stateFlds["potential_temperature"]);
  auto exnerLevelsHatView = make_view<double, 2>(hatFlds["exner_levels_minus_one"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto tHatView = make_view<double, 2>(hatFlds["air_temperature"]);

  idx_t lvls(hatFlds["air_temperature"].levels());
  idx_t lvlsm1 = lvls - 1;
  double exnerTopVal(0);
  double exnerTopHatVal(0);

  for (idx_t jn = 0; jn < tHatView.shape(0); ++jn) {
    // Passive code: Value above model top is assumed to be in hydrostatic balance.
    exnerTopVal = exnerLevelsView(jn, lvlsm1) -
      (constants::grav * (hlView(jn, lvls) - hlView(jn, lvlsm1))) /
      (constants::cp * thetaView(jn, lvlsm1));

    // Active code
    thetaHatView(jn, lvlsm1) += ( (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) * exnerTopVal +
      (hlView(jn, lvls)  - hView(jn, lvlsm1)) * exnerLevelsView(jn, lvlsm1) ) *
      tHatView(jn, lvlsm1) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));

    exnerTopHatVal = (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) *
      tHatView(jn, lvlsm1) * thetaView(jn, lvlsm1) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));

    exnerLevelsHatView(jn, lvlsm1) += (hlView(jn, lvls)  - hView(jn, lvlsm1)) *
      tHatView(jn, lvlsm1) * thetaView(jn, lvlsm1) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));

    tHatView(jn, lvlsm1) = 0.0;

    exnerLevelsHatView(jn, lvlsm1) += exnerTopHatVal;
    thetaHatView(jn, lvlsm1) += exnerTopHatVal * (exnerLevelsView(jn, lvlsm1) - exnerTopVal) /
        thetaView(jn, lvlsm1);
    exnerTopHatVal = 0.0;

    for (idx_t jl = lvls - 2; jl >= 0; --jl) {
      thetaHatView(jn, jl) += (
        (hView(jn, jl) - hlView(jn, jl)) * exnerLevelsView(jn, jl + 1) +
        (hlView(jn, jl + 1) - hView(jn, jl)) * exnerLevelsView(jn, jl) ) *
        tHatView(jn, jl) /
        (hlView(jn, jl + 1) - hlView(jn, jl));

      exnerLevelsHatView(jn, jl + 1) += (hView(jn, jl) - hlView(jn, jl)) *
        tHatView(jn, jl) * thetaView(jn, jl) /
        (hlView(jn, jl + 1) - hlView(jn, jl));

      exnerLevelsHatView(jn, jl) += (hlView(jn, jl + 1)  - hView(jn, jl)) *
        tHatView(jn, jl) * thetaView(jn, jl) /
        (hlView(jn, jl + 1) - hlView(jn, jl));

      tHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_air_temperature_ad()] ... exit" << std::endl;
}

}  // namespace mo
