/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <tuple>

#include "atlas/array.h"
#include "atlas/field/for_each.h"
#include "atlas/functionspace.h"

#include "mo/eval_relative_humidity.h"
#include "mo/eval_sat_vapour_pressure.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;
using atlas::util::Config;

namespace mo {

// --------------------------------------------------------------------------------------

bool eval_relative_humidity_nl(atlas::FieldSet & fields) {
  oops::Log::trace() << "[eval_relative_humidity_nl()] starting ..." << std::endl;

  bool cap_super_sat(false);

  if (fields["relative_humidity"].metadata().has("cap_super_sat")) {
    fields["relative_humidity"].metadata().get("cap_super_sat", cap_super_sat);
  }

  // Warning! Relative humidity in saber, vader and ufo has units percents (from 0 to 100),
  // while it should have units 1 (from 0 to 1) according to the CCPP convention.
  atlas::field::for_each_value(fields["specific_humidity"],
                               fields["qsat"],
                               fields["relative_humidity"],
                               [&](const double q, const double qsat, double& rh) {
      rh = fmax(q / qsat * 100.0, 0.0);
      rh = (cap_super_sat && (rh > 100)) ? 100.0 : rh;
  });

  oops::Log::trace() << "[eval_relative_humidity_nl()] ... exit" << std::endl;
  return true;
}

// --------------------------------------------------------------------------------------

void eval_relative_humidity_tl(atlas::FieldSet & incFlds,
                               const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_relative_humidity_tl()] starting ..." << std::endl;
  // * This assumes air pressure increments can be neglected (usual approximation)
  // * Only part of this is the tangent linear of eval_relative_humidity_nl.
  // * This neglects supersaturation in eval_relative_humidity_nl


  // State variables qsat and dlsvpdT can be calculated by eval_sat_vapour_pressure_nl
  // and evalSatSpecificHumidity from air_temperature and air_pressure inputs

  atlas::field::for_each_value(std::make_tuple(incFlds["specific_humidity"],
                                               incFlds["air_temperature"],
                                               incFlds["relative_humidity"],
                                               stateFlds["specific_humidity"],
                                               stateFlds["qsat"],
                                               stateFlds["dlsvpdT"]),
          [](const double qInc, const double tInc, double& rhInc,
             const double q, const double qsat, const double dlsvpdT) {
            rhInc = 100.0 * (qInc - q * dlsvpdT * tInc) / qsat;
          });

  oops::Log::trace() << "[eval_relative_humidity_tl()] ... exit" << std::endl;
}

// --------------------------------------------------------------------------------------

void eval_relative_humidity_ad(atlas::FieldSet & hatFlds,
                               const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_relative_humidity_ad()] starting ..." << std::endl;

  atlas::field::for_each_value(std::make_tuple(hatFlds["specific_humidity"],
                                               hatFlds["air_temperature"],
                                               hatFlds["relative_humidity"],
                                               stateFlds["specific_humidity"],
                                               stateFlds["qsat"],
                                               stateFlds["dlsvpdT"]),
          [](double& qHat, double& tHat, double& rhHat,
             const double q, const double qsat, const double dlsvpdT) {
            qHat += 100 * rhHat / qsat;
            tHat -= 100 * q * dlsvpdT * rhHat / qsat;
            rhHat = 0.0;
          });

  oops::Log::trace() << "[eval_relative_humidity_ad()] ... exit" << std::endl;
}

// --------------------------------------------------------------------------------------

void eval_relative_humidity_2m_tl(atlas::FieldSet & incFlds) {
  oops::Log::trace() << "[eval_relative_humidity_2m_tl()] starting ..." << std::endl;
  const auto rhIncView = make_view<double, 2>(incFlds["relative_humidity"]);
  auto rh2mIncView = make_view<double, 2>(incFlds["relative_humidity_2m"]);

  for (atlas::idx_t jn = 0; jn < rh2mIncView.shape(0); ++jn) {
      rh2mIncView(jn, 0) = rhIncView(jn, 0);
  }
  oops::Log::trace()<< "[eval_relative_humidity_2m_tl()] ... exit" << std::endl;
}

// --------------------------------------------------------------------------------------

void eval_relative_humidity_2m_ad(atlas::FieldSet & hatFlds) {
  oops::Log::trace() << "[eval_relative_humidity_2m_ad()] starting ..." << std::endl;
  auto rhHatView = make_view<double, 2>(hatFlds["relative_humidity"]);
  auto rh2mHatView = make_view<double, 2>(hatFlds["relative_humidity_2m"]);

  for (atlas::idx_t jn = 0; jn < rh2mHatView.shape(0); ++jn) {
    rhHatView(jn, 0) += rh2mHatView(jn, 0);
    rh2mHatView(jn, 0) = 0.0;
  }
  oops::Log::trace()<< "[eval_relative_humidity_2m_ad()] ... exit" << std::endl;
}

// --------------------------------------------------------------------------------------

}  // namespace mo
