/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

// TODO(Mayeul) Remove duplicate qtTemperature2qqclqcfTL in control2analysis_linearvarchange.h/cc
/// \details This is the moisture incrementing operator (MIO) that determines the
///          cloud ice water, cloud liquid water and specific humidity increments
///          from increments in total water and temperature
void eval_moisture_incrementing_operator_tl(atlas::FieldSet & incFlds,
                                            const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate qtTemperature2qqclqcfAD in control2analysis_linearvarchange.h/cc
/// \details This is the adjoint of the MIO
void eval_moisture_incrementing_operator_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate qqclqcf2qtTL in control2analysis_linearvarchange.h/cc
/// \details This is the inverse of the MIO, that determines the total water control
///          variable from specific humidity, cloud ice water and cloud liquid water
///          increments
void eval_total_water_tl(atlas::FieldSet & incFlds,
                         const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate qqclqcf2qtAD in control2analysis_linearvarchange.h/cc
/// \details This is the adjoint of eval_total_water_tl
void eval_total_water_ad(atlas::FieldSet & hatFlds,
                         const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate qqclqcf2qt in control2analysis_varchange.h/cc
/// \details This is the nonlinear version of eva_total_water_tl, that determines
///          the water vapor and cloud condensed water mixing ratios wrt moist
///          air and condensed water from the mixing ratios of each component
///          (water vapor, cloud ice and cloud liquid water).
void eval_total_water(atlas::FieldSet & stateFlds);

}  // namespace mo
