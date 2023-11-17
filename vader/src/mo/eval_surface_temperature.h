/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"

namespace mo {
/// \details This calculates the surface temperature from the atmospheric
///          temperature at the lowest level
void eval_surface_temperature_nl(atlas::FieldSet & stateFlds);

/// \details This calculates the surface temperature increment from the atmospheric
///          temperature increment at the lowest level
// TODO(Mayeul) Remove evalSurfaceTemperatureTL in model2geovals_linearvarchange.h/.cc
void eval_surface_temperature_tl(atlas::FieldSet & incFlds);

/// \details This calculates the adjoint of the atmospheric temperature
///          at the lowest level from the adjoint of the surface temperature
// TODO(Mayeul) Remove evalSurfaceTemperatureAD in model2geovals_linearvarchange.h/.cc
void eval_surface_temperature_adj(atlas::FieldSet & HatFlds);

}  // namespace mo
