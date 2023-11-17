/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"

#include "mo/constants.h"
#include "mo/model2geovals_linearvarchange.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;


namespace mo {

void evalAirPressureTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds) {
  // height at pressure levels
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  // height at theta levels
  const auto hView = make_view<const double, 2>(augStateFlds["height"]);

  const auto pIncView = make_view<double, 2>(incFlds["air_pressure_levels"]);
  auto pbarIncView = make_view<double, 2>(incFlds["air_pressure"]);

  double alpha_jl(0.0);

  atlas::idx_t lvls(incFlds["air_pressure"].levels());
  for (atlas::idx_t jn = 0; jn < pbarIncView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      alpha_jl = (hView(jn, jl) - hlView(jn, jl)) / (hlView(jn, jl+1) - hlView(jn, jl));
      pbarIncView(jn, jl) = (1.0 - alpha_jl) * pIncView(jn, jl) +
                            alpha_jl * pIncView(jn, jl+1);
    }
  }
}

void evalAirPressureAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds) {
  // height at pressure levels
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  // height at theta levels
  const auto hView = make_view<const double, 2>(augStateFlds["height"]);

  auto pbarHatView = make_view<double, 2>(hatFlds["air_pressure"]);
  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels"]);

  double alpha_jl(0.0);

  atlas::idx_t lvls(hatFlds["air_pressure"].levels());
  for (atlas::idx_t jn = 0; jn < pHatView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      alpha_jl = (hView(jn, jl) - hlView(jn, jl)) / (hlView(jn, jl+1) - hlView(jn, jl));
      pHatView(jn, jl) += (1.0 - alpha_jl) * pbarHatView(jn, jl);
      pHatView(jn, jl+1) += alpha_jl * pbarHatView(jn, jl);
      pbarHatView(jn, jl) = 0.0;
    }
  }
}

void evalRelativeHumidityTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds) {
  const auto qIncView = make_view<double, 2>(incFlds["specific_humidity"]);
  const auto temperIncView = make_view<double, 2>(incFlds["air_temperature"]);
  auto rhIncView = make_view<double, 2>(incFlds["relative_humidity"]);

  // qsat and dlnesdT are calculated by evalSatVaporPressure and evalSatSpecificHumidity
  // from air_temperature and air_pressure inputs
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto qsatView = make_view<const double, 2>(augStateFlds["qsat"]);
  const auto dlnesdT = make_view<const double, 2>(augStateFlds["dlsvpdT"]);

  atlas::idx_t lvls(incFlds["relative_humidity"].levels());
  for (atlas::idx_t jn = 0; jn < rhIncView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      rhIncView(jn, jl) = 100.0 * (qIncView(jn, jl) - qView(jn, jl) * dlnesdT(jn, jl)
                          * temperIncView(jn, jl)) / qsatView(jn, jl);
    }
  }
}

void evalRelativeHumidityAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds) {
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto temperHatView = make_view<double, 2>(hatFlds["air_temperature"]);
  auto rhHatView = make_view<double, 2>(hatFlds["relative_humidity"]);

  // qsat and dlnesdT are calculated by evalSatVaporPressure and evalSatSpecificHumidity
  // from air_temperature and air_pressure inputs
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto qsatView = make_view<const double, 2>(augStateFlds["qsat"]);
  const auto dlnesdT = make_view<const double, 2>(augStateFlds["dlsvpdT"]);

  atlas::idx_t lvls(hatFlds["relative_humidity"].levels());
  for (atlas::idx_t jn = 0; jn < rhHatView.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < lvls; ++jl) {
      qHatView(jn, jl) += 100.0 * rhHatView(jn, jl) / qsatView(jn, jl);
      temperHatView(jn, jl) -= 100.0 * qView(jn, jl) * dlnesdT(jn, jl) * rhHatView(jn, jl)
                            / qsatView(jn, jl);
      rhHatView(jn, jl) = 0.0;
    }
  }
}

void evalRelativeHumidity_2mTL(atlas::FieldSet & incFlds) {
  const auto rhIncView = make_view<double, 2>(incFlds["relative_humidity"]);
  auto rh2mIncView = make_view<double, 2>(incFlds["relative_humidity_2m"]);

  for (atlas::idx_t jn = 0; jn < rh2mIncView.shape(0); ++jn) {
      rh2mIncView(jn, 0) = rhIncView(jn, 0);
  }
}

void evalRelativeHumidity_2mAD(atlas::FieldSet & HatFlds) {
  auto rhHatView = make_view<double, 2>(HatFlds["relative_humidity"]);
  auto rh2mHatView = make_view<double, 2>(HatFlds["relative_humidity_2m"]);

  for (atlas::idx_t jn = 0; jn < rh2mHatView.shape(0); ++jn) {
    rhHatView(jn, 0) += rh2mHatView(jn, 0);
    rh2mHatView(jn, 0) = 0.0;
  }
}

void evalSurfacePressureTL(atlas::FieldSet & incFlds) {
  const auto pIncView = make_view<double, 2>(incFlds["air_pressure"]);
  auto pstarIncView = make_view<double, 2>(incFlds["surface_pressure"]);

  for (atlas::idx_t jn = 0; jn < pstarIncView.shape(0); ++jn) {
      pstarIncView(jn, 0) = pIncView(jn, 0);
  }
}

void evalSurfacePressureAD(atlas::FieldSet & HatFlds) {
  auto pHatView = make_view<double, 2>(HatFlds["air_pressure"]);
  auto pstarHatView = make_view<double, 2>(HatFlds["surface_pressure"]);

  for (atlas::idx_t jn = 0; jn < pstarHatView.shape(0); ++jn) {
    pHatView(jn, 0) += pstarHatView(jn, 0);
    pstarHatView(jn, 0) = 0.0;
  }
}

void evalSurfaceTemperatureTL(atlas::FieldSet & incFlds) {
  const auto temperIncView = make_view<double, 2>(incFlds["air_temperature"]);
  auto surfTemperIncView = make_view<double, 2>(incFlds["surface_temperature"]);

  for (atlas::idx_t jn = 0; jn < surfTemperIncView.shape(0); ++jn) {
      surfTemperIncView(jn, 0) = temperIncView(jn, 0);
  }
}

void evalSurfaceTemperatureAD(atlas::FieldSet & HatFlds) {
  auto temperHatView = make_view<double, 2>(HatFlds["air_temperature"]);
  auto surfTemperHatView = make_view<double, 2>(HatFlds["surface_temperature"]);

  for (atlas::idx_t jn = 0; jn < surfTemperHatView.shape(0); ++jn) {
    temperHatView(jn, 0) += surfTemperHatView(jn, 0);
    surfTemperHatView(jn, 0) = 0.0;
  }
}

}  // namespace mo
