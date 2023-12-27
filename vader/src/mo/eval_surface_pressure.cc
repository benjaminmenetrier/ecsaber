/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/eval_surface_pressure.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {

void eval_surface_pressure_tl(atlas::FieldSet & incFlds) {
  oops::Log::trace() << "[eval_surface_pressure_tl()] starting ..." << std::endl;
  const auto pIncView = make_view<double, 2>(incFlds["air_pressure_levels"]);
  auto pstarIncView = make_view<double, 2>(incFlds["surface_pressure"]);

  atlas_omp_parallel_for(atlas::idx_t jn = 0; jn < pstarIncView.shape(0); ++jn) {
      pstarIncView(jn, 0) = pIncView(jn, 0);
  }
  oops::Log::trace() << "[eval_surface_pressure_tl()] ... done" << std::endl;
}

void eval_surface_pressure_ad(atlas::FieldSet & HatFlds) {
  oops::Log::trace() << "[eval_surface_pressure_ad()] starting ..." << std::endl;
  auto pHatView = make_view<double, 2>(HatFlds["air_pressure_levels"]);
  auto pstarHatView = make_view<double, 2>(HatFlds["surface_pressure"]);

  atlas_omp_parallel_for(atlas::idx_t jn = 0; jn < pstarHatView.shape(0); ++jn) {
    pHatView(jn, 0) += pstarHatView(jn, 0);
    pstarHatView(jn, 0) = 0.0;
  }
  oops::Log::trace() << "[eval_surface_pressure_ad()] ... done" << std::endl;
}

}  // namespace mo
