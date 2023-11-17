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

/// \details This calculates the relative humidity increment from the air
///          temperature and humidity increment as:
///          rh' ~ (q'/qs - (q/qs) dlnesdT T')*100
///          This assumes air pressure increments can be neglected (usual approximation)
void evalRelativeHumidityTL(atlas::FieldSet & incFlds, const atlas::FieldSet & augStateFlds);

/// \details This calculates the adjoint of the air temperature and humidity
///          from the adjoint of the relative humidity
void evalRelativeHumidityAD(atlas::FieldSet & hatFlds, const atlas::FieldSet & augStateFlds);

/// \details This calculates the relative humidity increment at 2 m from the atmospheric
///          relative humidity increment at the lowest level
void evalRelativeHumidity_2mTL(atlas::FieldSet & incFlds);

/// \details This calculates the adjoint of the relative humidity
///          at the lowest level from the adjoint of the atmospheric relative humidity at 2 m
void evalRelativeHumidity_2mAD(atlas::FieldSet & HatFlds);

/// \details This calculates the surface pressure increment from the atmospheric
///          pressure increment at the lowest level
void evalSurfacePressureTL(atlas::FieldSet & incFlds);

/// \details This calculates the adjoint of the atmospheric pressure
///          at the lowest level from the adjoint of the surface pressure
void evalSurfacePressureAD(atlas::FieldSet & HatFlds);

/// \details This calculates the surface temperature increment from the atmospheric
///          temperature increment at the lowest level
void evalSurfaceTemperatureTL(atlas::FieldSet & incFlds);

/// \details This calculates the adjoint of the atmospheric temperature
///          at the lowest level from the adjoint of the surface temperature
void evalSurfaceTemperatureAD(atlas::FieldSet & HatFlds);

}  // namespace mo
