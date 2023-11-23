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
#include "atlas/util/Config.h"

#include "eckit/exception/Exceptions.h"

#include "mo/eval_ratio.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;
using atlas::util::Config;

namespace mo {

bool eval_ratio_to_second(atlas::FieldSet & fields, const std::vector<std::string> & vars)
{
  oops::Log::trace() << "[eval_ratio_to_second()] starting ..." << std::endl;

  if (vars.size() != 3) {
    throw eckit::Exception("Expected 3 variable names in vars", Here());
  }

  // vars[0] = m_x = [ mv | mci | mcl | m_r ]
  const auto ds_m_x  = make_view<const double, 2>(fields[vars[0]]);
  const auto ds_m_t  = make_view<const double, 2>(fields[vars[1]]);
  auto ds_tfield  = make_view<double, 2>(fields[vars[2]]);

  atlas_omp_parallel_for(idx_t jnode = 0; jnode < ds_tfield.shape(0); jnode++) {
    for (idx_t jlev = 0; jlev < ds_tfield.shape(1); jlev++) {
      ds_tfield(jnode, jlev) = ds_m_x(jnode, jlev) / ds_m_t(jnode, jlev);
    }
  }

  oops::Log::trace() << "[eval_ratio_to_second()] ... exit" << std::endl;
  return true;
}

// --------------------------------------------------------------------------------------

void eval_ratio_to_second_tl(atlas::FieldSet & incFields,
                             const atlas::FieldSet & fields,
                             const std::vector<std::string> & vars)
{
  oops::Log::trace() << "[eval_ratio_to_second_tl()] starting ..." << std::endl;

  if (vars.size() != 3) {
    throw eckit::Exception("Expected 3 variable names in vars", Here());
  }

  // vars[0] = m_x = [ mv | mci | mcl | m_r ]
  const auto inc_m_x  = make_view<const double, 2>(incFields[vars[0]]);
  const auto inc_m_t  = make_view<const double, 2>(incFields[vars[1]]);
  const auto m_x  = make_view<const double, 2>(fields[vars[0]]);
  const auto m_t  = make_view<const double, 2>(fields[vars[1]]);
  auto inc_ratio  = make_view<double, 2>(incFields[vars[2]]);

  atlas_omp_parallel_for(idx_t jnode = 0; jnode < m_t.shape(0); jnode++) {
    for (idx_t jlev = 0; jlev < m_t.shape(1); jlev++) {
      inc_ratio(jnode, jlev) = inc_m_x(jnode, jlev) / m_t(jnode, jlev)
                 - m_x(jnode, jlev) / (m_t(jnode, jlev) * m_t(jnode, jlev)) * inc_m_t(jnode, jlev);
    }
  }

  oops::Log::trace() << "[eval_ratio_to_second_tl()] ... exit" << std::endl;
}

// ----------------------------------------------------------------------------

void eval_ratio_to_second_ad(atlas::FieldSet & hatFields,
                             const atlas::FieldSet & fields,
                             const std::vector<std::string> & vars)
{
  oops::Log::trace() << "[eval_ratio_to_second_ad()] starting ..." << std::endl;

  if (vars.size() != 3) {
    throw eckit::Exception("Expected 3 variable names in vars", Here());
  }

  // vars[0] = m_x = [ mv | mci | mcl | m_r ]
  auto hat_m_x  = make_view<double, 2>(hatFields[vars[0]]);
  auto hat_m_t  = make_view<double, 2>(hatFields[vars[1]]);
  const auto m_x  = make_view<const double, 2>(fields[vars[0]]);
  const auto m_t  = make_view<const double, 2>(fields[vars[1]]);
  auto hat_ratio  = make_view<double, 2>(hatFields[vars[2]]);

  atlas_omp_parallel_for(idx_t jnode = 0; jnode < m_t.shape(0); jnode++) {
    for (idx_t jlev = 0; jlev < m_t.shape(1); jlev++) {
      hat_m_x(jnode, jlev) += hat_ratio(jnode, jlev) / m_t(jnode, jlev);
      hat_m_t(jnode, jlev) -= hat_ratio(jnode, jlev) * m_x(jnode, jlev)
                                / (m_t(jnode, jlev) * m_t(jnode, jlev));
      hat_ratio(jnode, jlev) = 0.0;
    }
  }

  oops::Log::trace() << "[eval_ratio_to_second_ad()] ... exit" << std::endl;
}

// --------------------------------------------------------------------------------------

}  // namespace mo
