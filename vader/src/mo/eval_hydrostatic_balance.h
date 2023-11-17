/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

// Almost a of duplicate hexner2PThetav in control2analysis_varchange.h/.cc
/// \details Hydrostatic balance from hydrostatic_exner_levels to virtual potential temperature
void eval_hydrobal_virtual_potential_temperature_nl(atlas::FieldSet & fields);

// TODO(MW) Delete duplicate hexner2ThetavTL in control2analysis_linearvarchange.h/.cc
/// \details Tangent linear approximation to the
///          transformation from hydrostatically-balanced exner (hydrostatic_exner_levels)
///          to virtual_potential_temperature
void eval_hydrobal_virtual_potential_temperature_tl(atlas::FieldSet & incFlds,
                                                    const atlas::FieldSet & augStateFlds);

// TODO(MW) Delete duplicate hexner2ThetavAD in control2analysis_linearvarchange.h/.cc
/// \details Adjoint of the tangent linear approximation to the
///          transformation from virtual potential temperature () to
///          hydrostatically-balanced exner (hexner)
void eval_hydrobal_virtual_potential_temperature_ad(atlas::FieldSet & hatFlds,
                                                    const atlas::FieldSet & augStateFlds);

/// \details Hydrostatic balance from virtual potential temperature and
///          air pressure levels minus one (at lowest level) to
///          hydrostatic_exner_levels
void eval_hydrobal_hydrostatic_exner_levels_nl(atlas::FieldSet & fields);

/// \details Tangent linear approximation to the transformation
///          from virtual potential temperature and
///          air pressure levels minus one (at lowest level) to
///          hydrostatic_exner_levels
void eval_hydrobal_hydrostatic_exner_levels_tl(atlas::FieldSet & incFlds,
                                               const atlas::FieldSet & augStateFlds);

/// \details Adjoint of the tangent linear approximation to the transformation
///          from virtual potential temperature and
///          air pressure levels minus one (at lowest level) to
///          hydrostatic_exner_levels
void eval_hydrobal_hydrostatic_exner_levels_ad(atlas::FieldSet & hatFlds,
                                               const atlas::FieldSet & augStateFlds);
}  // namespace mo
