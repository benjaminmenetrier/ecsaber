/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/eval_moisture_incrementing_operator.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace mo {

// ------------------------------------------------------------------------------------------------
// Note this is a copy of qtTemperature2qqclqcfTL
// TODO(Mayeul,MW) Remove duplicate in control2analysis_linearvarchange.cc/h
void eval_moisture_incrementing_operator_tl(atlas::FieldSet & incFlds,
                                            const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_moisture_incrementing_operator_tl()] starting ..."
                     << std::endl;
  auto qsatView = make_view<const double, 2>(augStateFlds["qsat"]);
  auto dlsvpdTView = make_view<const double, 2>(augStateFlds["dlsvpdT"]);
  auto cleffView = make_view<const double, 2>(augStateFlds["cleff"]);
  auto cfeffView = make_view<const double, 2>(augStateFlds["cfeff"]);
  auto qtIncView = make_view<const double, 2>(incFlds["qt"]);
  auto temperIncView = make_view<const double, 2>(incFlds["air_temperature"]);

  auto qclIncView = make_view<double, 2>
                    (incFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  auto qcfIncView = make_view<double, 2>
                    (incFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  auto qIncView = make_view<double, 2>(incFlds["specific_humidity"]);

  double maxCldInc;
  for (atlas::idx_t jn = 0; jn < incFlds["qt"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < incFlds["qt"].levels(); ++jl) {
        maxCldInc = qtIncView(jn, jl) - qsatView(jn, jl) *
            dlsvpdTView(jn, jl) * temperIncView(jn, jl);
        qclIncView(jn, jl) = cleffView(jn, jl) * maxCldInc;
        qcfIncView(jn, jl) = cfeffView(jn, jl) * maxCldInc;
        qIncView(jn, jl) = qtIncView(jn, jl) - qclIncView(jn, jl) - qcfIncView(jn, jl);
    }
  }
  oops::Log::trace() << "[eval_moisture_incrementing_operator_tl()] ... done"
                     << std::endl;
}

// ------------------------------------------------------------------------------------------------
// Note this is a copy of qtTemperature2qqclqcfAD
// TODO(Mayeul,MW) Remove duplicate in control2analysis_linearvarchange.cc/h
void eval_moisture_incrementing_operator_ad(atlas::FieldSet & hatFlds,
                             const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_moisture_incrementing_operator_ad()] starting ..."
                     << std::endl;
  auto qsatView = make_view<const double, 2>(augStateFlds["qsat"]);
  auto dlsvpdTView = make_view<const double, 2>(augStateFlds["dlsvpdT"]);
  auto cleffView = make_view<const double, 2>(augStateFlds["cleff"]);
  auto cfeffView = make_view<const double, 2>(augStateFlds["cfeff"]);

  auto temperHatView = make_view<double, 2>(hatFlds["air_temperature"]);
  auto qtHatView = make_view<double, 2>(hatFlds["qt"]);
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto qclHatView = make_view<double, 2>
                    (hatFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  auto qcfHatView = make_view<double, 2>
                    (hatFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);

  double qsatdlsvpdT;
  for (atlas::idx_t jn = 0; jn < hatFlds["qt"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < hatFlds["qt"].levels(); ++jl) {
      qsatdlsvpdT = qsatView(jn, jl) * dlsvpdTView(jn, jl);
      temperHatView(jn, jl) += ((cleffView(jn, jl) + cfeffView(jn, jl)) * qHatView(jn, jl)
                                - cleffView(jn, jl) * qclHatView(jn, jl)
                                - cfeffView(jn, jl) * qcfHatView(jn, jl)) * qsatdlsvpdT;
      qtHatView(jn, jl) += cleffView(jn, jl) * qclHatView(jn, jl)
              + cfeffView(jn, jl) * qcfHatView(jn, jl)
              + (1.0 - cleffView(jn, jl) - cfeffView(jn, jl))
              * qHatView(jn, jl);
      qHatView(jn, jl) = 0.0;
      qclHatView(jn, jl) = 0.0;
      qcfHatView(jn, jl) = 0.0;
    }
  }
  oops::Log::trace() << "[eval_moisture_incrementing_operator_ad()] ... done"
                     << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_total_water_tl(atlas::FieldSet & incFlds,
                         const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_total_water_tl()] starting ..." << std::endl;

  auto qIncView = make_view<const double, 2>(incFlds["specific_humidity"]);
  auto qclIncView = make_view<const double, 2>
                      (incFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  auto qcfIncView = make_view<const double, 2>
                      (incFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);

  auto qtIncView = make_view<double, 2>(incFlds["qt"]);

  for (atlas::idx_t jnode = 0; jnode < incFlds["qt"].shape(0); jnode++) {
    for (atlas::idx_t jlev = 0; jlev < incFlds["qt"].levels(); jlev++) {
      qtIncView(jnode, jlev) = qIncView(jnode, jlev)
                             + qclIncView(jnode, jlev)
                             + qcfIncView(jnode, jlev);
    }
  }
  oops::Log::trace() << "[eval_total_water_tl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_total_water_ad(atlas::FieldSet & hatFlds,
                         const atlas::FieldSet & augStateFlds) {
  oops::Log::trace() << "[eval_total_water_ad()] starting ..." << std::endl;
  auto qIncView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto qclIncView = make_view<double, 2>
                      (hatFlds["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  auto qcfIncView = make_view<double, 2>
                      (hatFlds["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  auto qtIncView = make_view<double, 2>(hatFlds["qt"]);

  for (atlas::idx_t jnode = 0; jnode < hatFlds["qt"].shape(0); jnode++) {
    for (atlas::idx_t jlev = 0; jlev < hatFlds["qt"].levels(); jlev++) {
      qIncView(jnode, jlev) += qtIncView(jnode, jlev);
      qclIncView(jnode, jlev) += qtIncView(jnode, jlev);
      qcfIncView(jnode, jlev) += qtIncView(jnode, jlev);
    }
  }

  qtIncView.assign(0.0);

  oops::Log::trace() << "[eval_total_water_ad()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_total_water(atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_total_water()] starting ..." << std::endl;
  eval_total_water_tl(stateFlds, stateFlds);
  oops::Log::trace() << "[eval_total_water()] ... done" << std::endl;
}


}  // namespace mo
