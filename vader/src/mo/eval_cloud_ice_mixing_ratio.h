/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

// TODO(Mayeul) Remove duplicate evalMassCloudIce in model2geovals_varchange.cc/.h
/// \brief function to evaluate the cloud ice water mixing ratio
///        wrt moist air and condensed water:
///   qci = m_ci/m_t
/// where ...
///   m_ci = cloud ice water mixing ratio wrt dry air
///   m_t  = moist air and condensed water mixing ratio wrt dry air
///
bool eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
        atlas::FieldSet & stateFlds);


void eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_tl(
        atlas::FieldSet & incFlds,
        const atlas::FieldSet & stateFlds);


void eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_ad(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds);

}  // namespace mo

