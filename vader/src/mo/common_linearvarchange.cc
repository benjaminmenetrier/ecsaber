/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/common_linearvarchange.h"
#include "mo/constants.h"
#include "mo/functions.h"

using atlas::array::make_view;
using atlas::idx_t;
using atlas::util::Config;

namespace mo {

/// \details Calculate the tangent linear of virtual potential temperature
///          from the specific humidity and the potential temperature.
void evalVirtualPotentialTemperatureTL(atlas::FieldSet & incFlds,
                                       const atlas::FieldSet & augStateFlds) {
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  const auto qIncView = make_view<const double, 2>(incFlds["specific_humidity"]);
  const auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  auto vthetaIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  auto fspace = incFlds["virtual_potential_temperature"].functionspace();

  auto evaluateVThetaTL = [&] (idx_t i, idx_t j) {
    vthetaIncView(i, j) = thetaView(i, j) * constants::c_virtual * qIncView(i, j) +
        thetaIncView(i, j) * (1.0 + constants::c_virtual * qView(i, j));
  };

  auto conf = Config("levels", incFlds["virtual_potential_temperature"].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateVThetaTL, conf);
}

/// \details Calculate the tangent linear of virtual potential temperature
///          from the specific humidity and the potential temperature.
void evalVirtualPotentialTemperatureAD(atlas::FieldSet & hatFlds,
                                       const atlas::FieldSet & augStateFlds) {
  const auto qView = make_view<const double, 2>(augStateFlds["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  auto qHatView = make_view<double, 2>(hatFlds["specific_humidity"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto vthetaHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);

  auto fspace = hatFlds["virtual_potential_temperature"].functionspace();

  auto evaluateVThetaAD = [&] (idx_t i, idx_t j) {
    qHatView(i, j) += thetaView(i, j) * constants::c_virtual * vthetaHatView(i, j);
    thetaHatView(i, j) +=  vthetaHatView(i, j) * (1.0 + constants::c_virtual * qView(i, j));
    vthetaHatView(i, j) = 0;
  };

  auto conf = Config("levels", hatFlds["virtual_potential_temperature"].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateVThetaAD, conf);
}


}  // namespace mo
