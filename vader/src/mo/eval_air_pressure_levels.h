/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

// TODO(Mayeul) Remove duplicate evalAirPressureLevels in common_varchange.cc/h
/// \details Vertical extrapolation of air pressure above model top
void eval_air_pressure_levels_nl(atlas::FieldSet & stateFlds);

void eval_air_pressure_levels_tl(atlas::FieldSet & incFlds, const atlas::FieldSet & stateFlds);

void eval_air_pressure_levels_ad(atlas::FieldSet & incFlds, const atlas::FieldSet & stateFlds);

}  // namespace mo
