/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"

namespace mo {

/// \details This calculates the surface pressure increment from the atmospheric
///          pressure increment at the lowest level
// TODO(Mayeul) Remove evalSurfacePressureTL in model2geovals_linearvarchange.h/.cc
void eval_surface_pressure_tl(atlas::FieldSet & incFlds);

/// \details This calculates the adjoint of the atmospheric pressure
///          at the lowest level from the adjoint of the surface pressure
// TODO(Mayeul) Remove evalSurfacePressureAD in model2geovals_linearvarchange.h/.cc
void eval_surface_pressure_ad(atlas::FieldSet & HatFlds);

}  // namespace mo
