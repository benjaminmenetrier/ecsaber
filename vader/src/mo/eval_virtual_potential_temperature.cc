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
#include "mo/eval_virtual_potential_temperature.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {

// ------------------------------------------------------------------------------------------------
void eval_virtual_potential_temperature_nl(atlas::FieldSet & fields) {
  oops::Log::trace() << "[eval_virtual_potential_temperature_nl()] starting ..." << std::endl;
  const auto qView = make_view<const double, 2>(fields["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(fields["potential_temperature"]);
  auto vThetaView = make_view<double, 2>(fields["virtual_potential_temperature"]);

  const atlas::idx_t n_levels = fields["virtual_potential_temperature"].levels();
  atlas_omp_parallel_for(atlas::idx_t ih = 0; ih < vThetaView.shape(0); ih++) {
    for (atlas::idx_t ilev = 0; ilev < n_levels; ilev++) {
      vThetaView(ih, ilev) = thetaView(ih, ilev)
              * (1.0 + constants::c_virtual * qView(ih, ilev));
    }
  }
  oops::Log::trace() << "[eval_virtual_potential_temperature_nl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_virtual_potential_temperature_tl(atlas::FieldSet & incFlds,
                                             const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_virtual_potential_temperature_tl()] starting ..." << std::endl;
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  const auto qIncView = make_view<const double, 2>(incFlds["specific_humidity"]);
  const auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  auto vThetaIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  const atlas::idx_t n_levels = incFlds["virtual_potential_temperature"].levels();
  atlas_omp_parallel_for(atlas::idx_t ih = 0; ih < vThetaIncView.shape(0); ih++) {
    for (atlas::idx_t ilev = 0; ilev < n_levels; ilev++) {
      vThetaIncView(ih, ilev) = thetaView(ih, ilev) * constants::c_virtual * qIncView(ih, ilev) +
                     thetaIncView(ih, ilev) * (1.0 + constants::c_virtual * qView(ih, ilev));
    }
  }
  oops::Log::trace() << "[eval_virtual_potential_temperature_tl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_virtual_potential_temperature_ad(atlas::FieldSet & hatFlds,
                                           const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_virtual_potential_temperature_ad()] starting ..." << std::endl;
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto vThetaHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);

  const atlas::idx_t n_levels = hatFlds["virtual_potential_temperature"].levels();
  atlas_omp_parallel_for(atlas::idx_t ih = 0; ih < vThetaHatView.shape(0); ih++) {
    for (atlas::idx_t ilev = 0; ilev < n_levels; ilev++) {
      qHatView(ih, ilev) += thetaView(ih, ilev) * constants::c_virtual * vThetaHatView(ih, ilev);
      thetaHatView(ih, ilev) +=  vThetaHatView(ih, ilev)
            * (1.0 + constants::c_virtual * qView(ih, ilev));
      vThetaHatView(ih, ilev) = 0;
    }
  }
  oops::Log::trace() << "[eval_virtual_potential_temperature_ad()] ... done" << std::endl;
}

}  // namespace mo
