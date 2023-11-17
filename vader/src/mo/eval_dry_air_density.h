/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

namespace mo {

/// \details Calculate the dry air density
void eval_dry_air_density_from_pressure_levels_minus_one_nl(atlas::FieldSet & stateFlds);

void eval_dry_air_density_from_pressure_levels_minus_one_tl(atlas::FieldSet & incFlds,
                                                      const atlas::FieldSet & stateFlds);

void eval_dry_air_density_from_pressure_levels_minus_one_ad(atlas::FieldSet & hatFlds,
                                                      const atlas::FieldSet & stateFlds);

void eval_dry_air_density_from_pressure_levels_nl(atlas::FieldSet & stateFlds);

void eval_dry_air_density_from_pressure_levels_tl(atlas::FieldSet & incFlds,
                                            const atlas::FieldSet & stateFlds);

void eval_dry_air_density_from_pressure_levels_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & stateFlds);

}  // namespace mo
