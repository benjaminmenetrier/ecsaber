/*
 * (C) Crown Copyright 2023 Met Office
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

// TODO(Mayeul) Remove duplicate evalSatVaporPressure in common_varchange.cc/h
/// \brief function to evaluate saturation water pressure (svp) [Pa]
/// the Atlas field in the argument must contain an inizialised air temperature field
/// and to have a defined svp field which is then calculated and returned as output
///
void eval_sat_vapour_pressure_nl(const std::string & svp_file, atlas::FieldSet & fields);

}
