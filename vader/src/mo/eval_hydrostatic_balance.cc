/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/constants.h"
#include "mo/eval_hydrostatic_balance.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

void eval_hydrobal_virtual_potential_temperature_nl(atlas::FieldSet & fields) {
  oops::Log::trace() <<
    "[eval_hydrobal_virtual_potential_temperature_nl] starting ..." << std::endl;
  auto hlView = make_view<const double, 2>(fields["height_levels"]);
  auto hexnerView = make_view<const double, 2>(fields["hydrostatic_exner_levels"]);
  auto thetavView = make_view<double, 2>(fields["virtual_potential_temperature"]);

  const atlas::idx_t n_levels(fields["virtual_potential_temperature"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < thetavView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      thetavView(ih, ilev) = - constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)) /
         (constants::cp * (hexnerView(ih, ilev+1) - hexnerView(ih, ilev)));
    }
  }
  oops::Log::trace() << "[eval_hydrobal_virtual_potential_temperature_nl] ... done" << std::endl;
}

void eval_hydrobal_virtual_potential_temperature_tl(atlas::FieldSet & incFlds,
                                             const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() <<
    "[eval_hydrobal_virtual_potential_temperature_tl()] starting ..." << std::endl;
  auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto hexnerIncView = make_view<const double, 2>(incFlds["hydrostatic_exner_levels"]);
  auto thetavIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  const idx_t n_levels = incFlds["virtual_potential_temperature"].levels();
  atlas_omp_parallel_for(idx_t ih = 0; ih < thetavIncView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      thetavIncView(ih, ilev) =
        (hexnerIncView(ih, ilev+1) - hexnerIncView(ih, ilev)) *
        (constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)) /
        (constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)));
    }
  }
  oops::Log::trace() << "[eval_hydrobal_virtual_potential_temperature_tl()] ... done" << std::endl;
}

void eval_hydrobal_virtual_potential_temperature_ad(atlas::FieldSet & hatFlds,
                                              const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() <<
    "[eval_hydrobal_virtual_potential_temperature_ad] starting ..." << std::endl;
  auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto hexnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  const idx_t n_levels_minus_one = hatFlds["virtual_potential_temperature"].levels() - 1;
  atlas_omp_parallel_for(idx_t ih = 0; ih < thetavHatView.shape(0); ++ih) {
    for (idx_t ilev = n_levels_minus_one; ilev > -1; --ilev) {
      hexnerHatView(ih, ilev+1) += thetavHatView(ih, ilev) *
        (constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)) /
        (constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)) );
      hexnerHatView(ih, ilev) -= thetavHatView(ih, ilev) *
        (constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)) /
        (constants::grav * (hlView(ih, ilev+1) - hlView(ih, ilev)));
      thetavHatView(ih, ilev) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_hydrobal_virtual_potential_temperature_ad] ... done" << std::endl;
}

// Note that this is a copy of evalHydrostaticExnerLevels
void eval_hydrobal_hydrostatic_exner_levels_nl(atlas::FieldSet & fields) {
  oops::Log::trace() << "[eval_hydrobal_hydrostatic_exner_levels_nl()] ... starting" << std::endl;
  auto rpView = make_view<const double, 2>(fields["height_levels"]);
  auto vThetaView = make_view<const double, 2>(fields["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(fields["air_pressure_levels"]);
  auto hexnerView = make_view<double, 2>(fields["hydrostatic_exner_levels"]);

  const atlas::idx_t n_levels(fields["hydrostatic_exner_levels"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < hexnerView.shape(0); ++ih) {
    hexnerView(ih, 0) = pow(pView(ih, 0) / constants::p_zero,
      constants::rd_over_cp);
    for (idx_t ilev = 1; ilev < n_levels; ++ilev) {
      hexnerView(ih, ilev) = hexnerView(ih, ilev-1) -
        (constants::grav * (rpView(ih, ilev) - rpView(ih, ilev-1))) /
        (constants::cp * vThetaView(ih, ilev-1));
    }
  }
  oops::Log::trace() << "[eval_hydrobal_hydrostatic_exner_levels_nl()] ... done" << std::endl;
}


void eval_hydrobal_hydrostatic_exner_levels_tl(atlas::FieldSet & incFlds,
                                               const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_hydrobal_hydrostatic_exner_levels_tl()] ... starting" << std::endl;

  auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(augStateFlds["air_pressure_levels"]);
  auto hexnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  auto thetavIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels"]);
  auto hexnerIncView = make_view<double, 2>(incFlds["hydrostatic_exner_levels"]);

  const idx_t n_levels_minus_one = incFlds["hydrostatic_exner_levels"].levels() - 1;
  atlas_omp_parallel_for(idx_t ih = 0; ih < hexnerIncView.shape(0); ++ih) {
    hexnerIncView(ih, 0) = constants::rd_over_cp *
      hexnerView(ih, 0) * pIncView(ih, 0) / pView(ih, 0);
    for (atlas::idx_t ilev = 0; ilev < n_levels_minus_one; ++ilev) {
      hexnerIncView(ih, ilev+1) = hexnerIncView(ih, ilev) +
        ((constants::grav * thetavIncView(ih, ilev) *
          (hlView(ih, ilev+1) - hlView(ih, ilev))) /
         (constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)));
    }
  }
  oops::Log::trace() << "[eval_hydrobal_hydrostatic_exner_levels_tl()] ... done" << std::endl;
}

void eval_hydrobal_hydrostatic_exner_levels_ad(atlas::FieldSet & hatFlds,
                                               const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_hydrobal_hydrostatic_exner_levels_ad()] ... starting" << std::endl;
  auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto pView = make_view<const double, 2>(augStateFlds["air_pressure_levels"]);
  auto hexnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);
  auto hexnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  const idx_t n_levels_minus_two = hatFlds["hydrostatic_exner_levels"].levels() - 2;
  atlas_omp_parallel_for(idx_t ih = 0; ih < hexnerHatView.shape(0); ++ih) {
    for (atlas::idx_t ilev = n_levels_minus_two; ilev >= 0; --ilev) {
      thetavHatView(ih, ilev) += ((constants::grav * hexnerHatView(ih, ilev+1) *
        (hlView(ih, ilev+1) - hlView(ih, ilev))) /
        (constants::cp * thetavView(ih, ilev) * thetavView(ih, ilev)));
      hexnerHatView(ih, ilev) += hexnerHatView(ih, ilev+1);
      hexnerHatView(ih, ilev+1) = 0.0;
    }
    pHatView(ih, 0) += constants::rd_over_cp *
      hexnerView(ih, 0) * hexnerHatView(ih, 0) / pView(ih, 0);
    hexnerHatView(ih, 0) = 0.0;
  }
  oops::Log::trace() << "[hydrobal_hydrostatic_exner_levels_ad()] ... done" << std::endl;
}
}  // namespace mo
