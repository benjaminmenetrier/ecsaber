/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

// TODO(Mayeul) Remove duplicate evalSpecificHumidity in model2geovals_varchange.cc/.h
/// \brief function to evaluate the water vapor mixing ratio wrt
///        moist air and condensed water (one possible definition
///        of specific humidity):
///   q = m_v/m_t
/// where ...
///   m_v  = water vapor mixing ratio wrt dry air
///   m_t  = moist air and condensed water mixing ratio wrt dry air
///
bool eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
        atlas::FieldSet & stateFlds);

void eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_tl(
        atlas::FieldSet & incFlds,
        const atlas::FieldSet & stateFlds);

void eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_ad(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds);


// TODO(Mayeul) Remove duplicate evalSpecificHumidityFromRH_2m in model2geovals_varchange.cc/.h
/// \brief function to evaluate the water vapor mixing ratio wrt moist air
///        and condensed water at 2m:
///   q2m = rh * qsat
/// where ...
///   q2m  = water vapor mixing ratio wrt moist air and condensed water at 2m above surface
///   rh   = relative humidity at 2m above surface
///   qsat = water vapor mixing ratio wrt moist air and condensed water assuming saturation
///
bool eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_at_2m_nl(
        atlas::FieldSet & stateFlds);

void eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_at_2m_tl(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds);

void eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_at_2m_ad(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds);

}  // namespace mo

