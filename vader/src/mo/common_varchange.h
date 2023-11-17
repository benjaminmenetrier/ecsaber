/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

namespace mo {

/// \brief function to evaluate saturation water pressure (svp) [Pa]
/// the Atlas field in the argument must contain an inizialised air temperature field
/// and to have a defined svp field which is then calculated and returned as output
///
bool evalSatVaporPressure(atlas::FieldSet & fields);

/// \brief function to evaluate saturation specific humidity (qsat)
/// Needs air pressure [Pa] and svp [Pa] Atlas fields and returns the qsat Atlas field
bool evalSatSpecificHumidity(atlas::FieldSet & fields);


/// \brief function to evaluate the 'air_pressure_levels' from
/// 'exner_levels_minus_one', 'potential_temperature', 'height_levels' fields
/// where ...
///   'air_pressure_levels' is the same as 'air_pressure_levels_minus_one'
/// except for the topmost level where the level is assumed to be in
/// hydrostatic balance (and capped to be >0).
///
bool evalAirPressureLevels(atlas::FieldSet & fields);


}  // namespace mo
