/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

namespace mo {

/// \details Tangent linear approximation to the
///          transformation from virtual potential temperature (thetav) to
///          hydrostatically-balanced exner (hydrostatic_exner_levels_minus_one)
void thetavP2HexnerTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds);

/// \details Adjoint of the tangent linear approximation to the
///          transformation from virtual potential temperature (thetav) to
///          hydrostatically-balanced exner (hexner)
void thetavP2HexnerAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds);

/// \details Tangent linear approximation to the
///          transformation from hydrostatically-balanced exner (hexner)
///          to virtual potential temperature (thetav)
void hexner2ThetavTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds);

/// \details Adjoint of the tangent linear approximation to the
///          transformation from virtual potential temperature (thetav) to
///          hydrostatically-balanced exner (hexner)
///          Note:: zeroing of hexnerHatView not done.
void hexner2ThetavAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds);

/// \details This calculates air temperature increments from "exner" and "theta"
///          with "exner" vertically interpolated to "theta" positions.
void evalAirTemperatureTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds);

/// \details This calculates the adjoint of the transform to air temperature increments.
void evalAirTemperatureAD(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds);

/// \details Calculate qT increment from the sum of q, qcl and qcf increments
void qqclqcf2qtTL(atlas::FieldSet & incFields, const atlas::FieldSet &);

/// \details This is the adjoint of qqclqcf2qtTL
void qqclqcf2qtAD(atlas::FieldSet & hatFields, const atlas::FieldSet &);

/// \details This is the moisture incrementing operator (MIO) that determines the cloud condensate
///          increments from increments in total water and temperature
void qtTemperature2qqclqcfTL(atlas::FieldSet & incFlds,
                             const atlas::FieldSet & augStateFlds);

/// \details This is the adjoint of the MIO
void qtTemperature2qqclqcfAD(atlas::FieldSet & hatFlds,
                             const atlas::FieldSet & augStateFlds);

/// \details This calculates the vertically-regressed geostrophic pressure increment field
///          in grid point space and adds it to the unbalanced pressure increment field to give
///          hydrostatic balance increments
void evalHydrostaticPressureTL(atlas::FieldSet & incFlds,
                               const atlas::FieldSet & augStateFlds);

/// \details This is the adjoint of the calculation of hydrostatic balance increment fields
void evalHydrostaticPressureAD(atlas::FieldSet & hatFlds,
                               const atlas::FieldSet & augStateFlds);

/// \details This calculates the hydrostatic exner field from the hydrostatic pressure
void evalHydrostaticExnerTL(atlas::FieldSet & incFlds,
                            const atlas::FieldSet & augStateFlds);

/// \details This is the adjoint of the calculation of hydrostatic exner increment fields
void evalHydrostaticExnerAD(atlas::FieldSet & hatFlds,
                            const atlas::FieldSet & augStateFlds);

/// \details This is calculates the moisture control variable and virtual potential temperature
///          increment fields
void evalMuThetavTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augState);

/// \details This is the adjoint of the calculation for the moisture control variable
///          and virtual potential temperature increments fields.
void evalMuThetavAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augState);

/// \details This calculates the total water (qt) and potential tempetature (theta)
///          increment fields.
void evalQtThetaTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augState);

/// \details This is the adjoint of the calculation to create the total water (qt)
///          and potential tempetature (theta) increment fields.
void evalQtThetaAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augState);

}  // namespace mo
