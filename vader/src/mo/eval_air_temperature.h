/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

// TODO(Mayeul): Remove duplicate evalAirTemperature in model2geovals_varchange.h/cc
void eval_air_temperature_nl(atlas::FieldSet & stateFlds);


// TODO(Mayeul): Remove duplicate evalAirTemperatureTL in control2analysis_linearvarchange.h/cc
void eval_air_temperature_tl(atlas::FieldSet & incFlds, const atlas::FieldSet & stateFlds);

// TODO(Mayeul): Remove duplicate evalAirTemperatureAD in control2analysis_linearvarchange.h/cc
void eval_air_temperature_ad(atlas::FieldSet & hatFlds, const atlas::FieldSet & stateFlds);

}  // namespace mo

