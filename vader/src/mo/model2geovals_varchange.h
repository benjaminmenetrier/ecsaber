/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"


/// \brief ++ Variable Tranforms Repository ++
///
/// \details This repository contains a set of functions associated with
/// 'transformations'.
/// The transformations are applied to Atlas data structures.
///
/// Here the term 'transformation' is used in a broad sense.
/// Examples of transformations are:
///   i) initialization of a field
///  ii) updating a field with an uniform value
/// iii) transformation applied to a variable associated with a given model (variable change)
///  iv) ...
///
/// This is only a temporary repository. In the long term, the Variable Transforms will be
/// stored and processed by the VAriable DErivation Repository (VADER) system.
///


namespace mo {

/// \brief function to initialize a field;
/// the field is stored into a Atlas data structure with rank=2
///
void initField_rank2(atlas::Field & field, const double value_init);


/// \brief function to update a field by using a single value,
/// the field is stored into a Atlas data structure with rank=2
///
void setUniformValue_rank2(atlas::Field & field, const double value);


/// \brief function to evaluate the 'total mass of moist air':
///   m_t = (1 + m_v + m_ci + m_cl + m_r)
/// where ...
///   m_t  = total mass of moist air
///   m_v  = mixing ratio of water vapour
///   m_ci = mixing ratio of cloud ice
///   m_cl = mixing ratio of cloud liquid
///   m_r  = mixing ratio of rain
///
bool evalTotalMassMoistAir(atlas::FieldSet & fields);

/// \brief function to evaluate the quantity:
///   qx = m_x/m_t
/// where ...
///   m_x = [ mv | mci | mcl | m_r ]
///   m_t  = total mass of moist air
///
bool evalRatioToMt(atlas::FieldSet & fields);


/// \brief function to evaluate the 'specific humidity':
///   q = m_v/m_t
/// where ...
///   m_v  = mixing ratio of water vapour
///   m_t  = total mass of moist air
///
bool evalSpecificHumidity(atlas::FieldSet & fields);

/// \brief function to evaluate the 'total relative humidity':
///   rh = (q+qcl+qci+qrain)/qsat*100
/// where ...
///   q  = specific humidity
///   qcl = specific cloud water
///   qci = specific cloud ice
///   qrain = specific rain
///   qsat  = saturated specific humidity
///
bool evalTotalRelativeHumidity(atlas::FieldSet & fields);

/// \brief function to evaluate the 'mass content of cloud ice' in atmosphere layer:
///   qci = m_ci/m_t
/// where ...
///   m_ci = mixing ratio of cloud ice
///   m_t  = total mass of moist air
///
bool evalMassCloudIce(atlas::FieldSet & fields);


/// \brief function to evaluate the 'mass content of liquid water' in atmosphere layer:
///   qcl = m_cl/m_t
/// where ...
///   m_cl = mixing ratio of cloud liquid
///   m_t  = total mass of moist air
///
bool evalMassCloudLiquid(atlas::FieldSet & fields);


/// \brief function to evaluate the 'mass content of rain' in atmosphere layer:
///   qrain = m_r/m_t
/// where ...
///   m_r = mixing ratio of rain
///   m_t = total mass of moist air
///
bool evalMassRain(atlas::FieldSet & fields);


/// \brief function to evaluate the 'air temperature':
///   air_temperature = theta x exner
/// where ...
///   theta = theta, air potential temperature (LFRic names)
///   exner = exner_in_wth, exner pressure on wtheta points (LFRic names)
///
/// note that ...
/// 'theta' and 'exner' are on the same levels
///
bool evalAirTemperature(atlas::FieldSet & fields);

/// \brief function to evaluate the 'specific humidity at two meters above surface':
///   q2m = rh * qsat
/// where ...
///   q2m  = specific humidity at two meters above surface
///   rh   = relative humidity at two meters above surface
///   qsat = saturation specific humidity
///
bool evalSpecificHumidityFromRH_2m(atlas::FieldSet & fields);


/// \brief function to evaluate 'surf_param_a' and 'surf_param_b';
/// 'surf_param_a' and 'surf_param_b' are  two 'derived' fields used to evaluate
/// the background pressure at the observation height
///
///          surf_param_a = f(height_levels[0], t_msh)
///          surf_param_b = f(t_msh, pressure_level_minus_one[0])
/// where ...
///     height_levels[0] = model surface height
///            t_msh = temperature at model surface height
///     pressure_level_minus_one[0] = pressure at model surface height
///
/// input data (fields) needed for the calculations:
///     height = height on theta model levels (not including surface)
///        (it includes metadata "boundary_layer_index" )
///     height_levels = height on rho model levels
///     pressure_levels_minus_one = pressure on rho model levels
///     specific_humidity
///
bool evalParamAParamB(atlas::FieldSet & fields);

/// \brief function to calculate Geopotential height field
///        on rho levels as equal to the height field on rho levels.
///        This is a good approximation in the lower atmosphere
///        where the acceleration due to gravity (g) is approximatey
///        equal to g0, the globally averaged accelaration due
///        to gravity at the Earth's surface
///
bool evalGpHeightMinusOne(atlas::FieldSet & fields);

}  // namespace mo
