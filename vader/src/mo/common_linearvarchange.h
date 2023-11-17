/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"

namespace mo {

/// \details Calculate the virtual potential temperature
///          from the specific humidity and the potential temperature.
void evalVirtualPotentialTemperatureTL(atlas::FieldSet & incFlds,
                                       const atlas::FieldSet & augStateFlds);

/// \details This calculates the adjoint of the transform to air temperature increments.
void evalVirtualPotentialTemperatureAD(atlas::FieldSet & incFlds,
                                       const atlas::FieldSet & augStateFlds);

}  // namespace mo
