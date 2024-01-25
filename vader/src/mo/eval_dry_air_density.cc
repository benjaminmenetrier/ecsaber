/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"

#include "mo/constants.h"
#include "mo/eval_dry_air_density.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

// -------------------------------------------------------------------------------------------------

/// \details Calculate the dry air density from virtual potential temperature
///          and air pressure (using air_pressure_levels_minus_one)
void eval_dry_air_density_from_pressure_levels_minus_one_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_nl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto rhoView = make_view<double, 2>(stateFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = stateFlds["dry_air_density_levels_minus_one"].levels();

  double h_minus_hm1;
  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < stateFlds["dry_air_density_levels_minus_one"].shape(0); ++jn) {
    rhoView(jn, 0) = pView(jn, 0) *
      pow(pView(jn, 0) / constants::p_zero, - constants::rd_over_cp) /
      (constants::rd * vptView(jn, 0));

    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hm1 = hView(jn, jl) - hView(jn, jl-1);
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoView(jn, jl) = pView(jn, jl) *
        pow(pView(jn, jl) / constants::p_zero, - constants::rd_over_cp)
        * h_minus_hm1 / (constants::rd * (
        h_minus_hl * vptView(jn, jl-1) + hl_minus_hm1 * vptView(jn, jl)));
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_nl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_minus_one_tl(atlas::FieldSet & incFlds,
                                                      const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_tl()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels_minus_one"]);
  auto vptIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto rhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = incFlds["dry_air_density_levels_minus_one"].levels();

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoIncView.shape(0); ++jn) {
    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoIncView(jn, jl) = rhoView(jn, jl) * ( (1.0 - constants::rd_over_cp) *
        pIncView(jn, jl) / pView(jn, jl) -
        ((hl_minus_hm1 * vptIncView(jn, jl) + h_minus_hl * vptIncView(jn, jl-1)) /
         (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1))));
    }

    rhoIncView(jn, 0) = rhoView(jn, 0) * (
      (1.0 - constants::rd_over_cp) * pIncView(jn, 0) / pView(jn, 0) -
      vptIncView(jn, 0)/ vptView(jn, 0));
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_minus_one_ad(atlas::FieldSet & hatFlds,
                                                     const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_ad()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels_minus_one"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels_minus_one"]);
  auto vptHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto rhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoHatView.shape(0); ++jn) {
    pHatView(jn, 0) += rhoView(jn, 0) * rhoHatView(jn, 0) * (1.0 - constants::rd_over_cp) /
      pView(jn, 0);
    vptHatView(jn, 0) -= rhoView(jn, 0) * rhoHatView(jn, 0) /
      vptView(jn, 0);
    rhoHatView(jn, 0) = 0.0;

    for (idx_t jl = hatFlds["dry_air_density_levels_minus_one"].levels()-1; jl >= 1; --jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      pHatView(jn, jl) += rhoView(jn, jl) * rhoHatView(jn, jl) * (1.0 - constants::rd_over_cp) /
        pView(jn, jl);
      vptHatView(jn, jl) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        hl_minus_hm1 /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      vptHatView(jn, jl-1) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        h_minus_hl /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      rhoHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_ad()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

/// \details Calculate the dry air density from virtual potential temperature
///          and air pressure (using air_pressure_levels)
void eval_dry_air_density_from_pressure_levels_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_nl()] starting ..."
                     << std::endl;
  const auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(stateFlds["height"]);
  const auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto rhoView = make_view<double, 2>(stateFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = stateFlds["dry_air_density_levels_minus_one"].levels();

  double h_minus_hm1;
  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < stateFlds["dry_air_density_levels_minus_one"].shape(0); ++jn) {
    rhoView(jn, 0) = pView(jn, 0) *
      pow(pView(jn, 0) / constants::p_zero, - constants::rd_over_cp) /
      (constants::rd * vptView(jn, 0));

    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hm1 = hView(jn, jl) - hView(jn, jl-1);
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoView(jn, jl) = pView(jn, jl) *
        pow(pView(jn, jl) / constants::p_zero, - constants::rd_over_cp)
        * h_minus_hm1 / (constants::rd * (
        h_minus_hl * vptView(jn, jl-1) + hl_minus_hm1 * vptView(jn, jl)));
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_nl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_tl(atlas::FieldSet & incFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_tl()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels"]);
  auto vptIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto rhoIncView = make_view<double, 2>(incFlds["dry_air_density_levels_minus_one"]);
  const idx_t numLevels = incFlds["dry_air_density_levels_minus_one"].levels();

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoIncView.shape(0); ++jn) {
    for (idx_t jl = 1; jl < numLevels; ++jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      rhoIncView(jn, jl) = rhoView(jn, jl) * ( (1.0 - constants::rd_over_cp) *
        pIncView(jn, jl) / pView(jn, jl) -
        ((hl_minus_hm1 * vptIncView(jn, jl) + h_minus_hl * vptIncView(jn, jl-1)) /
         (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1))));
    }

    rhoIncView(jn, 0) = rhoView(jn, 0) * (
      (1.0 - constants::rd_over_cp) * pIncView(jn, 0) / pView(jn, 0) -
      vptIncView(jn, 0)/ vptView(jn, 0));
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_pressure_levels_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_ad()] starting ..."
                     << std::endl;
  auto hlView = make_view<const double, 2>(stateFlds["height_levels"]);
  auto hView = make_view<const double, 2>(stateFlds["height"]);
  auto vptView = make_view<const double, 2>(stateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(stateFlds["air_pressure_levels"]);
  auto rhoView = make_view<const double, 2>(stateFlds["dry_air_density_levels_minus_one"]);

  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);
  auto vptHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto rhoHatView = make_view<double, 2>(hatFlds["dry_air_density_levels_minus_one"]);

  double h_minus_hl;
  double hl_minus_hm1;

  for (idx_t jn = 0; jn < rhoHatView.shape(0); ++jn) {
    pHatView(jn, 0) += rhoView(jn, 0) * rhoHatView(jn, 0) * (1.0 - constants::rd_over_cp) /
      pView(jn, 0);
    vptHatView(jn, 0) -= rhoView(jn, 0) * rhoHatView(jn, 0) /
      vptView(jn, 0);
    rhoHatView(jn, 0) = 0.0;

    for (idx_t jl = hatFlds["dry_air_density_levels_minus_one"].levels()-1; jl >= 1; --jl) {
      h_minus_hl = hView(jn, jl) - hlView(jn, jl);
      hl_minus_hm1 = hlView(jn, jl) - hView(jn, jl-1);
      pHatView(jn, jl) += rhoView(jn, jl) * rhoHatView(jn, jl) * (1.0 - constants::rd_over_cp) /
        pView(jn, jl);
      vptHatView(jn, jl) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        hl_minus_hm1 /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      vptHatView(jn, jl-1) -= rhoView(jn, jl) * rhoHatView(jn, jl) *
        h_minus_hl /
        (hl_minus_hm1 * vptView(jn, jl) + h_minus_hl * vptView(jn, jl-1));
      rhoHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_ad()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace mo
