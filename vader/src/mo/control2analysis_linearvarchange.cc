/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "mo/constants.h"
#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"

#include "atlas/array/MakeView.h"

using atlas::array::make_view;

namespace mo {

void thetavP2HexnerTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds) {
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  const auto thetavView = make_view<const double, 2>(
    augStateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(augStateFlds["air_pressure_levels_minus_one"]);
  const auto hexnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  const auto thetavIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  const auto pIncView = make_view<const double, 2>(incFlds["air_pressure_levels_minus_one"]);
  auto hexnerIncView = make_view<double, 2>(incFlds["hydrostatic_exner_levels"]);

  for (atlas::idx_t jn = 0; jn < incFlds["hydrostatic_exner_levels"].shape(0); ++jn) {
    hexnerIncView(jn, 0) = constants::rd_over_cp *
      hexnerView(jn, 0) * pIncView(jn, 0) / pView(jn, 0);

    for (atlas::idx_t jl = 1; jl < incFlds["hydrostatic_exner_levels"].levels(); ++jl) {
      hexnerIncView(jn, jl) = hexnerIncView(jn, jl-1) +
        ((constants::grav * thetavIncView(jn, jl-1) *
          (hlView(jn, jl) - hlView(jn, jl-1))) /
         (constants::cp * thetavView(jn, jl-1) * thetavView(jn, jl-1)));
    }
  }
}

void thetavP2HexnerAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds) {
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  const auto thetavView = make_view<const double, 2>(
    augStateFlds["virtual_potential_temperature"]);
  const auto pView = make_view<const double, 2>(augStateFlds["air_pressure_levels_minus_one"]);
  const auto hexnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto pHatView = make_view<double, 2>(hatFlds["air_pressure_levels_minus_one"]);
  auto hexnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  for (atlas::idx_t jn = 0; jn < hatFlds["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = hatFlds["hydrostatic_exner_levels"].levels() - 1; jl > 0; --jl) {
      thetavHatView(jn, jl-1) = thetavHatView(jn, jl-1) +
        ((constants::grav * hexnerHatView(jn, jl) *
        (hlView(jn, jl) - hlView(jn, jl-1))) /
        (constants::cp * thetavView(jn, jl-1) * thetavView(jn, jl-1)));

      hexnerHatView(jn, jl-1) = hexnerHatView(jn, jl-1) +
        hexnerHatView(jn, jl);
      hexnerHatView(jn, jl) = 0.0;
    }
    pHatView(jn, 0) = pHatView(jn, 0) +
      constants::rd_over_cp *
      hexnerView(jn, 0) * hexnerHatView(jn, 0) / pView(jn, 0);
    hexnerHatView(jn, 0) = 0.0;
  }
}

void hexner2ThetavTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds) {
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  const auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  const auto hexnerIncView = make_view<const double, 2>(incFlds["hydrostatic_exner_levels"]);
  auto thetavIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  atlas::idx_t levels = incFlds["virtual_potential_temperature"].levels();
  for (atlas::idx_t jn = 0; jn < incFlds["virtual_potential_temperature"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < levels; ++jl) {
      thetavIncView(jn, jl) =
        (hexnerIncView(jn, jl+1) - hexnerIncView(jn, jl)) *
        (constants::cp * thetavView(jn, jl) * thetavView(jn, jl)) /
        (constants::grav * (hlView(jn, jl+1) - hlView(jn, jl)));
    }
  }
}

void hexner2ThetavAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds) {
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  const auto thetavView = make_view<const double, 2>(augStateFlds["virtual_potential_temperature"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto hexnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  atlas::idx_t levelsm1 = hatFlds["virtual_potential_temperature"].levels()-1;
  for (atlas::idx_t jn = 0; jn < hatFlds["virtual_potential_temperature"].shape(0); ++jn) {
    for (atlas::idx_t jl = levelsm1; jl > -1; --jl) {
      hexnerHatView(jn, jl+1) += thetavHatView(jn, jl) *
        (constants::cp * thetavView(jn, jl) * thetavView(jn, jl)) /
        (constants::grav * (hlView(jn, jl+1) - hlView(jn, jl)) );
      hexnerHatView(jn, jl) -= thetavHatView(jn, jl) *
        (constants::cp * thetavView(jn, jl) * thetavView(jn, jl)) /
        (constants::grav * (hlView(jn, jl+1) - hlView(jn, jl)));
      thetavHatView(jn, jl) = 0.0;
    }
  }
}


/// \details This calculates air temperature increments.
void evalAirTemperatureTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds) {
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(augStateFlds["height"]);
  const auto exnerLevelsView = make_view<const double, 2>(augStateFlds["exner_levels_minus_one"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  const auto exnerLevelsIncView = make_view<const double, 2>(incFlds["exner_levels_minus_one"]);
  const auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  auto tIncView = make_view<double, 2>(incFlds["air_temperature"]);

  atlas::idx_t lvls(incFlds["air_temperature"].levels());
  atlas::idx_t lvlsm1 = lvls - 1;
  double exnerTopVal;
  double exnerTopIncVal;

  // Active code
  for (atlas::idx_t jn = 0; jn < tIncView.shape(0); ++jn) {
    // Passive code: Value above model top is assumed to be in hydrostatic balance.
    exnerTopVal = exnerLevelsView(jn, lvlsm1) -
      (constants::grav * (hlView(jn, lvls) - hlView(jn, lvlsm1))) /
      (constants::cp * thetaView(jn, lvlsm1));

    // Active code;
    for (atlas::idx_t jl= 0; jl < lvls - 1; ++jl) {
       tIncView(jn, jl) = (
         ( (hView(jn, jl) - hlView(jn, jl)) * exnerLevelsView(jn, jl + 1) +
           (hlView(jn, jl+1)  - hView(jn, jl)) * exnerLevelsView(jn, jl) ) *
           thetaIncView(jn, jl) +
         ( (hView(jn, jl) - hlView(jn, jl)) * exnerLevelsIncView(jn, jl + 1) +
           (hlView(jn, jl+1)  - hView(jn, jl)) * exnerLevelsIncView(jn, jl) ) *
         thetaView(jn, jl) ) /
         (hlView(jn, jl+1) - hlView(jn, jl));
    }

    exnerTopIncVal = exnerLevelsIncView(jn, lvlsm1) +
      thetaIncView(jn, lvlsm1) * (exnerLevelsView(jn, lvlsm1) - exnerTopVal) /
      thetaView(jn, lvlsm1);

    tIncView(jn, lvlsm1) = (
      ( (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) * exnerTopVal +
        (hlView(jn, lvls)  - hView(jn, lvlsm1)) * exnerLevelsView(jn, lvlsm1) ) *
         thetaIncView(jn, lvlsm1) +
      ( (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) * exnerTopIncVal +
        (hlView(jn, lvls)  - hView(jn, lvlsm1)) * exnerLevelsIncView(jn, lvlsm1) ) *
         thetaView(jn, lvlsm1)) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));
  }
}


/// \details This calculates air temperature increments.
void evalAirTemperatureAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds) {
  const auto hlView = make_view<const double, 2>(augStateFlds["height_levels"]);
  const auto hView = make_view<const double, 2>(augStateFlds["height"]);
  const auto exnerLevelsView = make_view<const double, 2>(augStateFlds["exner_levels_minus_one"]);
  const auto thetaView = make_view<const double, 2>(augStateFlds["potential_temperature"]);
  auto exnerLevelsHatView = make_view<double, 2>(hatFlds["exner_levels_minus_one"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto tHatView = make_view<double, 2>(hatFlds["air_temperature"]);

  atlas::idx_t lvls(hatFlds["air_temperature"].levels());
  atlas::idx_t lvlsm1 = lvls - 1;
  double exnerTopVal(0);
  double exnerTopHatVal(0);

  for (atlas::idx_t jn = 0; jn < tHatView.shape(0); ++jn) {
    // Passive code: Value above model top is assumed to be in hydrostatic balance.
    exnerTopVal = exnerLevelsView(jn, lvlsm1) -
      (constants::grav * (hlView(jn, lvls) - hlView(jn, lvlsm1))) /
      (constants::cp * thetaView(jn, lvlsm1));

    // Active code
    thetaHatView(jn, lvlsm1) += ( (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) * exnerTopVal +
      (hlView(jn, lvls)  - hView(jn, lvlsm1)) * exnerLevelsView(jn, lvlsm1) ) *
      tHatView(jn, lvlsm1) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));

    exnerTopHatVal = (hView(jn, lvlsm1) - hlView(jn, lvlsm1)) *
      tHatView(jn, lvlsm1) * thetaView(jn, lvlsm1) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));

    exnerLevelsHatView(jn, lvlsm1) += (hlView(jn, lvls)  - hView(jn, lvlsm1)) *
      tHatView(jn, lvlsm1) * thetaView(jn, lvlsm1) /
      (hlView(jn, lvls) - hlView(jn, lvlsm1));

    tHatView(jn, lvlsm1) = 0.0;

    exnerLevelsHatView(jn, lvlsm1) += exnerTopHatVal;
    thetaHatView(jn, lvlsm1) += exnerTopHatVal * (exnerLevelsView(jn, lvlsm1) - exnerTopVal) /
        thetaView(jn, lvlsm1);
    exnerTopHatVal = 0.0;

    for (atlas::idx_t jl = lvls - 2; jl >= 0; --jl) {
      thetaHatView(jn, jl) += (
        (hView(jn, jl) - hlView(jn, jl)) * exnerLevelsView(jn, jl + 1) +
        (hlView(jn, jl + 1) - hView(jn, jl)) * exnerLevelsView(jn, jl) ) *
        tHatView(jn, jl) /
        (hlView(jn, jl + 1) - hlView(jn, jl));

      exnerLevelsHatView(jn, jl + 1) += (hView(jn, jl) - hlView(jn, jl)) *
        tHatView(jn, jl) * thetaView(jn, jl) /
        (hlView(jn, jl + 1) - hlView(jn, jl));

      exnerLevelsHatView(jn, jl) += (hlView(jn, jl + 1)  - hView(jn, jl)) *
        tHatView(jn, jl) * thetaView(jn, jl) /
        (hlView(jn, jl + 1) - hlView(jn, jl));

      tHatView(jn, jl) = 0.0;
    }
  }
}


void qqclqcf2qtTL(atlas::FieldSet & incFields, const atlas::FieldSet &) {
  qqclqcf2qt(incFields);
}

void qqclqcf2qtAD(atlas::FieldSet & hatFields, const atlas::FieldSet &) {
  auto qHatView = make_view<double, 2>(hatFields["specific_humidity"]);
  auto qclHatView = make_view<double, 2>
                    (hatFields["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  auto qcfHatView = make_view<double, 2>
                    (hatFields["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  auto qtHatView = make_view<double, 2>(hatFields["qt"]);

  for (atlas::idx_t jn = 0; jn < hatFields["qt"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < hatFields["qt"].levels(); ++jl) {
      qHatView(jn, jl) += qtHatView(jn, jl);
      qclHatView(jn, jl) += qtHatView(jn, jl);
      qcfHatView(jn, jl) += qtHatView(jn, jl);
      qtHatView(jn, jl) = 0.0;
    }
  }
}

void qtTemperature2qqclqcfTL(atlas::FieldSet & incFlds,
                             const atlas::FieldSet & augStateFlds) {
  const auto qsatView = make_view<const double, 2>(augStateFlds["qsat"]);
  const auto dlsvpdTView = make_view<const double, 2>(augStateFlds["dlsvpdT"]);
  const auto cleffView = make_view<const double, 2>(augStateFlds["cleff"]);
  const auto cfeffView = make_view<const double, 2>(augStateFlds["cfeff"]);

  const auto qtIncView = make_view<const double, 2>(incFlds["qt"]);
  const auto temperIncView = make_view<const double, 2>(incFlds["air_temperature"]);
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
}

void qtTemperature2qqclqcfAD(atlas::FieldSet & hatFlds,
                             const atlas::FieldSet & augStateFlds) {
  const auto qsatView = make_view<const double, 2>(augStateFlds["qsat"]);
  const auto dlsvpdTView = make_view<const double, 2>(augStateFlds["dlsvpdT"]);
  const auto cleffView = make_view<const double, 2>(augStateFlds["cleff"]);
  const auto cfeffView = make_view<const double, 2>(augStateFlds["cfeff"]);

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
}


void evalHydrostaticPressureTL(atlas::FieldSet & incFlds,
                               const atlas::FieldSet & augStateFlds) {
  const auto gPIncView = make_view<const double, 2>(
    incFlds["geostrophic_pressure_levels_minus_one"]);
  const auto uPIncView = make_view<const double, 2>(
    incFlds["unbalanced_pressure_levels_minus_one"]);

  const auto pView = make_view<const double, 2>(augStateFlds["air_pressure_levels"]);
  // First index of interpWeightView is horizontal index, the second is bin index here

  atlas::idx_t levels = incFlds["geostrophic_pressure_levels_minus_one"].levels();
  auto hPIncView = make_view<double, 2>(incFlds["hydrostatic_pressure_levels"]);
  for (atlas::idx_t jn = 0; jn < incFlds["hydrostatic_pressure_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < levels; ++jl) {
      hPIncView(jn, jl) = uPIncView(jn, jl);
    }
  }

  if (augStateFlds.has("interpolation_weights")) {
    // Bins Vertical regression matrix stored in one field
    // B = (vertical regression matrix bin_0)
    //     (vertical regression matrix bin_1)
    //     (          ...                   )
    //     (vertical regression matrix bin_m)
    // Since each matrix is square we can easily infer the bin index from the row index
    // First index of vertRegView is bin_index * number of levels + level index,
    // the second is number of levels associated with matrix column.
    const auto vertRegView = make_view<const double, 2>(
      augStateFlds["vertical_regression_matrices"]);
    const auto interpWeightView = make_view<const double, 2>(
      augStateFlds["interpolation_weights"]);
    atlas::idx_t nBins = augStateFlds["interpolation_weights"].shape(1);

    for (atlas::idx_t jn = 0; jn < incFlds["hydrostatic_pressure_levels"].shape(0); ++jn) {
      for (atlas::idx_t jl = 0; jl < levels; ++jl) {
        for (atlas::idx_t b = 0; b < nBins; ++b) {
          for (atlas::idx_t jl2 = 0; jl2 < levels; ++jl2) {
            hPIncView(jn, jl) += interpWeightView(jn, b) *
                                 vertRegView(b * levels + jl, jl2) *
                                 gPIncView(jn, jl2);
          }
        }
      }
    }
  }

  for (atlas::idx_t jn = 0; jn < incFlds["hydrostatic_pressure_levels"].shape(0); ++jn) {
    hPIncView(jn, levels) =
      hPIncView(jn, levels-1) *
      std::pow(pView(jn, levels-1) / pView(jn, levels), constants::rd_over_cp - 1.0);
  }
}


void evalHydrostaticPressureAD(atlas::FieldSet & hatFlds,
                               const atlas::FieldSet & augStateFlds) {
  auto gPHatView = make_view<double, 2>(hatFlds["geostrophic_pressure_levels_minus_one"]);
  auto uPHatView = make_view<double, 2>(hatFlds["unbalanced_pressure_levels_minus_one"]);

  const auto pView = make_view<const double, 2>(augStateFlds["air_pressure_levels"]);
  auto hPHatView = make_view<double, 2>(hatFlds["hydrostatic_pressure_levels"]);

  atlas::idx_t levels = hatFlds["geostrophic_pressure_levels_minus_one"].levels();
  for (atlas::idx_t jn = 0; jn < hatFlds["hydrostatic_pressure_levels"].shape(0); ++jn) {
    hPHatView(jn, levels - 1) +=
     hPHatView(jn, levels) *
     std::pow(pView(jn, levels-1) / pView(jn, levels), constants::rd_over_cp - 1.0);
    hPHatView(jn, levels) = 0.0;
  }

  if (augStateFlds.has("interpolation_weights")) {
    // Bins Vertical regression matrix stored in one field
    // B = (vertical regression matrix bin_0)
    //     (vertical regression matrix bin_1)
    //     (          ...                   )
    //     (vertical regression matrix bin_m)
    // Since each matrix is square we can easily infer the bin index from the row index
    // First index of vertRegView is bin_index * number of levels + level index,
    //     the second is level index
    const auto vertRegView = make_view<const double, 2>(
      augStateFlds["vertical_regression_matrices"]);
    const auto interpWeightView = make_view<const double, 2>(
      augStateFlds["interpolation_weights"]);
    atlas::idx_t nBins = augStateFlds["interpolation_weights"].shape(1);

    for (atlas::idx_t jn = 0; jn < hatFlds["hydrostatic_pressure_levels"].shape(0); ++jn) {
      for (atlas::idx_t jl = levels - 1; jl >= 0; --jl) {
        for (atlas::idx_t b = nBins -1; b >= 0; --b) {
          for (atlas::idx_t jl2 = levels - 1; jl2 >= 0; --jl2) {
            gPHatView(jn, jl2) += interpWeightView(jn, b) *
                                  vertRegView(b * levels + jl, jl2) *
                                  hPHatView(jn, jl);
          }
        }
        uPHatView(jn, jl) += hPHatView(jn, jl);
        hPHatView(jn, jl) = 0.0;
      }
    }
  }

  for (atlas::idx_t jn = 0; jn < hatFlds["hydrostatic_pressure_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = levels - 1; jl >= 0; --jl) {
      uPHatView(jn, jl) += hPHatView(jn, jl);
      hPHatView(jn, jl) = 0.0;
    }
  }
}

/// \details This calculates the hydrostatic exner field from the hydrostatic pressure
void evalHydrostaticExnerTL(atlas::FieldSet & incFlds,
                            const atlas::FieldSet & augStateFlds) {
  const auto pView = make_view<const double, 2>(augStateFlds["hydrostatic_pressure_levels"]);
  const auto exnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  const auto pIncView = make_view<const double, 2>(incFlds["hydrostatic_pressure_levels"]);
  auto exnerIncView = make_view<double, 2>(incFlds["hydrostatic_exner_levels"]);

  atlas::idx_t levels = incFlds["hydrostatic_exner_levels"].levels();
  for (atlas::idx_t jn = 0; jn < incFlds["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < levels; ++jl) {
      exnerIncView(jn, jl) = pIncView(jn, jl) *
        (constants::rd_over_cp * exnerView(jn, jl)) /
        pView(jn, jl);
    }
  }
}

/// \details This is the adjoint of the calculation of hydrostatic exner increments
void evalHydrostaticExnerAD(atlas::FieldSet & hatFlds,
                            const atlas::FieldSet & augStateFlds) {
  const auto pView = make_view<const double, 2>(augStateFlds["hydrostatic_pressure_levels"]);
  const auto exnerView = make_view<const double, 2>(augStateFlds["hydrostatic_exner_levels"]);
  auto pHatView = make_view<double, 2>(hatFlds["hydrostatic_pressure_levels"]);
  auto exnerHatView = make_view<double, 2>(hatFlds["hydrostatic_exner_levels"]);

  atlas::idx_t levels = hatFlds["hydrostatic_exner_levels"].levels();
  for (atlas::idx_t jn = 0; jn < hatFlds["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < levels; ++jl) {
      pHatView(jn, jl) += exnerHatView(jn, jl) *
        (constants::rd_over_cp * exnerView(jn, jl)) /
        pView(jn, jl);
      exnerHatView(jn, jl) = 0.0;
    }
  }
}


/// \details This is function calculates the linear moisture
///          control variable (muInc and thetavInc) from (thetaInc) and (qtInc)
///          We are ignoring the scaled pressure contribution to mu, because we have
///          found in the past that it gives no benefit and that its contribution
///          is small.
void evalMuThetavTL(atlas::FieldSet & incFlds,  const atlas::FieldSet & augState) {
  const auto muRow1Column1View = make_view<const double, 2>(augState["muRow1Column1"]);
  const auto muRow1Column2View = make_view<const double, 2>(augState["muRow1Column2"]);
  const auto muRow2Column1View = make_view<const double, 2>(augState["muRow2Column1"]);
  const auto muRow2Column2View = make_view<const double, 2>(augState["muRow2Column2"]);
  const auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  const auto qtIncView = make_view<const double, 2>(incFlds["qt"]);
  auto muIncView = make_view<double, 2>(incFlds["mu"]);
  auto thetavIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  for (atlas::idx_t jn = 0; jn < incFlds["mu"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < incFlds["mu"].levels(); ++jl) {
      muIncView(jn, jl) = muRow1Column1View(jn, jl)  * qtIncView(jn, jl)
                        + muRow1Column2View(jn, jl)  * thetaIncView(jn, jl);
      thetavIncView(jn, jl) = muRow2Column1View(jn, jl)  * qtIncView(jn, jl)
                            + muRow2Column2View(jn, jl)  * thetaIncView(jn, jl);
    }
  }
}


void evalMuThetavAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augState) {
  const auto muRow1Column1View = make_view<const double, 2>(augState["muRow1Column1"]);
  const auto muRow1Column2View = make_view<const double, 2>(augState["muRow1Column2"]);
  const auto muRow2Column1View = make_view<const double, 2>(augState["muRow2Column1"]);
  const auto muRow2Column2View = make_view<const double, 2>(augState["muRow2Column2"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto qtHatView = make_view<double, 2>(hatFlds["qt"]);
  auto muHatView = make_view<double, 2>(hatFlds["mu"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);

  for (atlas::idx_t jn = 0; jn < hatFlds["mu"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < hatFlds["mu"].levels(); ++jl) {
      thetaHatView(jn, jl) += muRow2Column2View(jn, jl) * thetavHatView(jn, jl);
      qtHatView(jn, jl) += muRow2Column1View(jn, jl) * thetavHatView(jn, jl);
      thetaHatView(jn, jl) += muRow1Column2View(jn, jl) * muHatView(jn, jl);
      qtHatView(jn, jl) += muRow1Column1View(jn, jl) * muHatView(jn, jl);
      thetavHatView(jn, jl) = 0.0;
      muHatView(jn, jl) = 0.0;
    }
  }
}


void evalQtThetaTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augState) {
  // Using Cramer's rule to calculate inverse.
  const auto muRecipDeterView = make_view<const double, 2>(augState["muRecipDeterminant"]);
  const auto muRow1Column1View = make_view<const double, 2>(augState["muRow1Column1"]);
  const auto muRow1Column2View = make_view<const double, 2>(augState["muRow1Column2"]);
  const auto muRow2Column1View = make_view<const double, 2>(augState["muRow2Column1"]);
  const auto muRow2Column2View  = make_view<const double, 2>(augState["muRow2Column2"]);
  const auto muIncView = make_view<const double, 2>(incFlds["mu"]);
  const auto thetavIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto qtIncView = make_view<double, 2>(incFlds["qt"]);
  auto thetaIncView = make_view<double, 2>(incFlds["potential_temperature"]);

  for (atlas::idx_t jn = 0; jn < incFlds["mu"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < incFlds["mu"].levels(); ++jl) {
      // VAR equivalent in Var_UpPFtheta_qT.f90 for thetaIncView
      // (beta2 * muA * theta_v' +   beta1 * mu') /
      // (alpha1 * beta2 * muA - alpha2 * muA * beta1)
      thetaIncView(jn, jl) =  muRecipDeterView(jn, jl) * (
                             muRow1Column1View(jn, jl) * thetavIncView(jn, jl)
                           - muRow2Column1View(jn, jl) * muIncView(jn, jl) );

      // VAR equivalent in Var_UpPFtheta_qT.f90 for qtIncView
      // (alpha1 * mu_v' -   alpha2 * muA * thetav') /
      // (alpha1 * beta2 * muA - alpha2 * muA * beta1)
      qtIncView(jn, jl) =  muRecipDeterView(jn, jl) * (
                           muRow2Column2View(jn, jl) * muIncView(jn, jl) -
                           muRow1Column2View(jn, jl) * thetavIncView(jn, jl) );
    }
  }
}


void evalQtThetaAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augState) {
  const auto muRecipDeterView = make_view<const double, 2>(augState["muRecipDeterminant"]);
  const auto muRow1Column1View = make_view<const double, 2>(augState["muRow1Column1"]);
  const auto muRow1Column2View = make_view<const double, 2>(augState["muRow1Column2"]);
  const auto muRow2Column1View = make_view<const double, 2>(augState["muRow2Column1"]);
  const auto muRow2Column2View  = make_view<const double, 2>(augState["muRow2Column2"]);
  auto qtHatView = make_view<double, 2>(hatFlds["qt"]);
  auto muHatView = make_view<double, 2>(hatFlds["mu"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);

  for (atlas::idx_t jn = 0; jn < hatFlds["mu"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < hatFlds["mu"].levels(); ++jl) {
      thetavHatView(jn, jl) += muRecipDeterView(jn, jl) *
                               muRow1Column1View(jn, jl) * thetaHatView(jn, jl);
      muHatView(jn, jl) -= muRecipDeterView(jn, jl) *
                           muRow2Column1View(jn, jl) * thetaHatView(jn, jl);
      thetavHatView(jn, jl) -= muRecipDeterView(jn, jl) *
                               muRow1Column2View(jn, jl) * qtHatView(jn, jl);
      muHatView(jn, jl) += muRecipDeterView(jn, jl) *
                           muRow2Column2View(jn, jl) * qtHatView(jn, jl);
      thetaHatView(jn, jl) = 0.0;
      qtHatView(jn, jl) = 0.0;
    }
  }
}

}  // namespace mo
