/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

/// \details This calculates the air pressure increment fields on theta (full) levels
///          from vertical interpolation of air pressure increment fields
///          on pressure (half) levels:
///          p'(i) = ((h(i) - h(i-1/2))p'(i+1/2) + (h(i+1/2) - h(i))p'(i-1/2))
///                  / (h(i+1/2) - h(i-1/2))
void evalAirPressureTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds);

/// \details This is the adjoint of the calculation to create the air pressure increment fields
///          on theta (full) levels from vertical interpolation of air pressure increment fields
///          on pressure (half) levels.
void evalAirPressureAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds);

}  // namespace mo
