/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/eval_surface_temperature.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {

void eval_surface_temperature_nl(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_surface_temperature_nl()] starting ..." << std::endl;
  const auto temperView = make_view<double, 2>(stateFlds["air_temperature"]);
  auto surfTemperView = make_view<double, 2>(stateFlds["surface_temperature"]);

  atlas_omp_parallel_for(atlas::idx_t jn = 0; jn < surfTemperView.shape(0); ++jn) {
    surfTemperView(jn, 0) = temperView(jn, 0);
  }
  oops::Log::trace() << "[eval_surface_temperature_nl()] ... done" << std::endl;
}

void eval_surface_temperature_tl(atlas::FieldSet & incFlds) {
  oops::Log::trace() << "[eval_surface_temperature_tl()] starting ..." << std::endl;
  const auto temperIncView = make_view<double, 2>(incFlds["air_temperature"]);
  auto surfTemperIncView = make_view<double, 2>(incFlds["surface_temperature"]);

  atlas_omp_parallel_for(atlas::idx_t jn = 0; jn < surfTemperIncView.shape(0); ++jn) {
      surfTemperIncView(jn, 0) = temperIncView(jn, 0);
  }
  oops::Log::trace() << "[eval_surface_temperature_tl()] ... done" << std::endl;
}

void eval_surface_temperature_adj(atlas::FieldSet & HatFlds) {
  oops::Log::trace() << "[eval_surface_temperature_adj()] starting ..." << std::endl;
  auto temperHatView = make_view<double, 2>(HatFlds["air_temperature"]);
  auto surfTemperHatView = make_view<double, 2>(HatFlds["surface_temperature"]);

  atlas_omp_parallel_for(atlas::idx_t jn = 0; jn < surfTemperHatView.shape(0); ++jn) {
    temperHatView(jn, 0) += surfTemperHatView(jn, 0);
    surfTemperHatView(jn, 0) = 0.0;
  }
  oops::Log::trace() << "[eval_surface_temperature_adj()] ... done" << std::endl;
}

}  // namespace mo
