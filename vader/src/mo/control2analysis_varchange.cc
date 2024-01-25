/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/Config.h"
#include "mo/constants.h"
#include "mo/control2analysis_varchange.h"
#include "mo/functions.h"

using atlas::array::make_view;
using atlas::util::Config;
using atlas::idx_t;

namespace mo {


void hexner2PThetav(atlas::FieldSet & fields) {
  const auto rpView = make_view<const double, 2>(fields["height_levels"]);
  const auto hexnerView = make_view<const double, 2>(fields["hydrostatic_exner_levels"]);
  auto pView = make_view<double, 2>(fields["air_pressure_levels_minus_one"]);
  auto vthetaView = make_view<double, 2>(fields["virtual_potential_temperature"]);
  const idx_t numLevels = fields["hydrostatic_exner_levels"].levels();

  for (idx_t jn = 0; jn < fields["hydrostatic_exner_levels"].shape(0); ++jn) {
    pView(jn, 0) = constants::p_zero * pow(hexnerView(jn, 0), (constants::cp / constants::rd));

    for (idx_t jl = 1; jl < numLevels; ++jl) {
      vthetaView(jn, jl) = -constants::grav * (rpView(jn, jl) - rpView(jn, jl-1)) /
         (constants::cp * (hexnerView(jn, jl) - hexnerView(jn, jl-1)));
    }
    vthetaView(jn, 0) = vthetaView(jn, 1);
  }
}

void evalVirtualPotentialTemperature(atlas::FieldSet & fields) {
  const auto qView = make_view<const double, 2>(fields["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(fields["potential_temperature"]);
  auto vthetaView = make_view<double, 2>(fields["virtual_potential_temperature"]);

  auto fspace = fields["virtual_potential_temperature"].functionspace();

  auto evaluateVTheta = [&] (idx_t i, idx_t j) {
    vthetaView(i, j) = thetaView(i, j) * (1.0 + constants::c_virtual * qView(i, j)); };

  auto conf = Config("levels", fields["virtual_potential_temperature"].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateVTheta, conf);
}

/// \details Calculate the hydrostatic exner pressure (on levels)
///          using air_pressure_minus_one and virtual potential temperature.
void evalHydrostaticExnerLevels(atlas::FieldSet & fields) {
  const auto rpView = make_view<const double, 2>(fields["height_levels"]);
  const auto vthetaView = make_view<const double, 2>(fields["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(fields["air_pressure_levels_minus_one"]);
  auto hexnerView = make_view<double, 2>(fields["hydrostatic_exner_levels"]);
  const idx_t numLevels = fields["hydrostatic_exner_levels"].levels();

  for (idx_t jn = 0; jn < fields["hydrostatic_exner_levels"].shape(0); ++jn) {
    hexnerView(jn, 0) = pow(pView(jn, 0) / constants::p_zero,
      constants::rd_over_cp);
    for (idx_t jl = 1; jl < numLevels; ++jl) {
      hexnerView(jn, jl) = hexnerView(jn, jl-1) -
        (constants::grav * (rpView(jn, jl) - rpView(jn, jl-1))) /
        (constants::cp * vthetaView(jn, jl-1));
    }
  }
}


/// \details Calculate the hydrostatic pressure (on levels)
///           from hydrostatic exner
void evalHydrostaticPressureLevels(atlas::FieldSet & fields) {
  const auto hexnerView = make_view<double, 2>(fields["hydrostatic_exner_levels"]);
  auto hpView = make_view<double, 2>(fields["hydrostatic_pressure_levels"]);
  const idx_t numLevels = fields["hydrostatic_pressure_levels"].levels();

  for (idx_t jn = 0; jn < fields["hydrostatic_pressure_levels"].shape(0); ++jn) {
    for (idx_t jl = 0; jl < numLevels; ++jl) {
       hpView(jn, jl) = constants::p_zero *
         pow(hexnerView(jn, jl), 1.0 / constants::rd_over_cp);
    }
  }
}


/// \details Calculate qT increment from the sum of q, qcl and qcf increments
void qqclqcf2qt(atlas::FieldSet & fields) {
  const auto qIncView = make_view<const double, 2>(fields["specific_humidity"]);
  const auto qclIncView = make_view<const double, 2>
                    (fields["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  const auto qcfIncView = make_view<const double, 2>
                    (fields["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  auto qtIncView = make_view<double, 2>(fields["qt"]);
  const idx_t numLevels = fields["specific_humidity"].levels();

  for (atlas::idx_t jn = 0; jn < fields["specific_humidity"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < numLevels; ++jl) {
      qtIncView(jn, jl) = qIncView(jn, jl) + qclIncView(jn, jl) + qcfIncView(jn, jl);
    }
  }
}

/// \details Calculate exner pressure levels
///          from air_pressure_levels_minus_one and using hydrostatic balance relation
///          for topmost level
void evalExnerPressureLevels(atlas::FieldSet & fields) {
  const auto exnerMinusOneView = make_view<const double, 2>(fields["exner_levels_minus_one"]);
  // Note that it is unclear whether this should be virtual_potential_temperature
  // or potential_temperature in this case. Either way the difference will be tiny since
  // the amount of moisture at a model top is tiny.
  const auto vthetaView = make_view<const double, 2>(fields["virtual_potential_temperature"]);
  const auto hlView = make_view<const double, 2>(fields["height_levels"]);
  auto exnerView = make_view<double, 2>(fields["exner_pressure_levels"]);

  idx_t levels(fields["exner_pressure_levels"].levels());
  for (idx_t jn = 0; jn < fields["exner_pressure_levels"].shape(0); ++jn) {
    for (idx_t jl = 1; jl < levels - 1; ++jl) {
      exnerView(jn, jl) = exnerMinusOneView(jn, jl);
    }

    exnerView(jn, levels - 1) = exnerView(jn, levels - 2) -
      (constants::grav * (hlView(jn, levels - 1) - hlView(jn, levels - 2))) /
      (constants::cp * vthetaView(jn, levels - 2));

    exnerView(jn, levels - 1) = exnerView(jn, levels-1) > 0.0 ?
      exnerView(jn, levels - 1) : constants::deps;
  }
}


void evalMoistureControlDependencies(atlas::FieldSet & fields) {
  const auto qtView = make_view<const double, 2>(fields["qt"]);
  const auto qView = make_view<const double, 2>(fields["specific_humidity"]);
  const auto thetaView = make_view<const double, 2>(fields["potential_temperature"]);
  const auto exnerView = make_view<const double, 2>(fields["exner"]);
  const auto dlsvpdTView = make_view<const double, 2>(fields["dlsvpdT"]);
  const auto qsatView = make_view<const double, 2>(fields["qsat"]);
  const auto muAView = make_view<const double, 2>(fields["muA"]);
  const auto muH1View = make_view<const double, 2>(fields["muH1"]);

  // this is effectively the (2x2) matrix = A
  //  (mu')       = A (qt')     where A is
  //  (theta_v')      (theta')
  //
  //  ( muA/qsat    - (muA/qsat) muH1 qT exner_bar dlsvpdT )
  //  (                                                 )
  //  (c_v theta q     (1+ cv) q                        )
  auto muRow1Column1View = make_view<double, 2>(fields["muRow1Column1"]);
  auto muRow1Column2View = make_view<double, 2>(fields["muRow1Column2"]);
  auto muRow2Column1View = make_view<double, 2>(fields["muRow2Column1"]);
  auto muRow2Column2View = make_view<double, 2>(fields["muRow2Column2"]);
  auto muRecipDeterminantView = make_view<double, 2>(fields["muRecipDeterminant"]);

  const idx_t numLevels = fields["potential_temperature"].levels();

  // the comments below are there to allow checking with the VAR code.
  for (atlas::idx_t jn = 0; jn < fields["potential_temperature"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < numLevels; ++jl) {
      muRow1Column1View(jn, jl) = muAView(jn, jl) / qsatView(jn, jl);  // beta2 * muA
      muRow1Column2View(jn, jl) = -  qtView(jn, jl)  * muH1View(jn, jl)
        * exnerView(jn, jl) * dlsvpdTView(jn, jl) * muRow1Column1View(jn, jl);
      // alpha2 * muA
      muRow2Column1View(jn, jl) = constants::c_virtual * thetaView(jn, jl);   // beta1
      muRow2Column2View(jn, jl) = 1.0 + constants::c_virtual * qView(jn, jl);  // alpha1
      muRecipDeterminantView(jn, jl) = 1.0 /(
        muRow2Column2View(jn, jl) * muRow1Column1View(jn, jl)
        - muRow1Column2View(jn, jl) * muRow2Column1View(jn, jl));
           // 1/( alpha1 * beta2 * muA - alpha2 * muA * beta1)
    }
  }
}

}  // namespace mo
