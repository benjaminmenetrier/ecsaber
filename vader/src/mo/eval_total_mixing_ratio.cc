/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/eval_ratio.h"
#include "mo/eval_total_mixing_ratio.h"

#include "oops/util/Logger.h"

namespace mo {

using atlas::array::make_view;
using atlas::idx_t;

// --------------------------------------------------------------------------------------

bool eval_total_mixing_ratio_nl(
        atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_nl()] starting ..."
          << std::endl;

  const auto ds_m_v  = make_view<const double, 2>(stateFlds["m_v"]);
  const auto ds_m_ci = make_view<const double, 2>(stateFlds["m_ci"]);
  const auto ds_m_cl = make_view<const double, 2>(stateFlds["m_cl"]);
  const auto ds_m_r  = make_view<const double, 2>(stateFlds["m_r"]);
  auto ds_m_t  = make_view<double, 2>(stateFlds["m_t"]);

  atlas_omp_parallel_for(idx_t jnode = 0; jnode < ds_m_t.shape(0); jnode++) {
    for (idx_t jlev = 0; jlev < ds_m_t.shape(1); jlev++) {
      ds_m_t(jnode, jlev) = 1.0 + ds_m_v(jnode, jlev) + ds_m_ci(jnode, jlev) +
                            ds_m_cl(jnode, jlev) + ds_m_r(jnode, jlev);
    }
  }

  oops::Log::trace()
          << "[eval_total_mixing_ratio_nl()] ... exit"
          << std::endl;
  return true;
}

// --------------------------------------------------------------------------------------

void eval_total_mixing_ratio_tl(
        atlas::FieldSet & incFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_tl()] starting ..."
          << std::endl;

  const auto inc_m_v  = make_view<const double, 2>(incFlds["m_v"]);
  const auto inc_m_ci = make_view<const double, 2>(incFlds["m_ci"]);
  const auto inc_m_cl = make_view<const double, 2>(incFlds["m_cl"]);
  const auto inc_m_r  = make_view<const double, 2>(incFlds["m_r"]);
  auto inc_m_t  = make_view<double, 2>(incFlds["m_t"]);

  atlas_omp_parallel_for(idx_t jnode = 0; jnode < inc_m_t.shape(0); jnode++) {
    for (idx_t jlev = 0; jlev < inc_m_t.shape(1); jlev++) {
      inc_m_t(jnode, jlev) = inc_m_v(jnode, jlev) + inc_m_ci(jnode, jlev) +
                             inc_m_cl(jnode, jlev) + inc_m_r(jnode, jlev);
    }
  }

  oops::Log::trace()
          << "[eval_total_mixing_ratio_tl()] ... exit"
          << std::endl;
}
// --------------------------------------------------------------------------------------

void eval_total_mixing_ratio_ad(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_total_mixing_ratio_ad()] starting ..."
          << std::endl;

  auto hat_m_v  = make_view<double, 2>(hatFlds["m_v"]);
  auto hat_m_ci = make_view<double, 2>(hatFlds["m_ci"]);
  auto hat_m_cl = make_view<double, 2>(hatFlds["m_cl"]);
  auto hat_m_r  = make_view<double, 2>(hatFlds["m_r"]);
  auto hat_m_t  = make_view<double, 2>(hatFlds["m_t"]);

  atlas_omp_parallel_for(idx_t jnode = 0; jnode < hat_m_t.shape(0); jnode++) {
    for (idx_t jlev = 0; jlev < hat_m_t.shape(1); jlev++) {
      hat_m_v(jnode, jlev) += hat_m_t(jnode, jlev);
      hat_m_ci(jnode, jlev) += hat_m_t(jnode, jlev);
      hat_m_cl(jnode, jlev) += hat_m_t(jnode, jlev);
      hat_m_r(jnode, jlev) += hat_m_t(jnode, jlev);
      hat_m_t(jnode, jlev) = 0.0;
    }
  }

  oops::Log::trace()
          << "[eval_total_mixing_ratio_ad()] ... exit"
          << std::endl;
}

// --------------------------------------------------------------------------------------
}  // namespace mo
