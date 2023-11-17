/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once


#include "atlas/field/FieldSet.h"

namespace mo {

// TODO(Mayeul) Remove duplicate evalMoistureControlDependencies in control2analysis_varchange.h/cc
/// \details The linear moisture control variable trajectory fields
///          All fields needs to be allocated beforehand.
///          Inputs: qt
///                  specific humidity
///                  potential temperature
///                  exner
///                  dlsvpdT
///                  qsat
///                  muA (normalisation factor based on statistics)
///                  muH1 (statistical factor that affects whether moisture
///                        control variable is acting as qt
///                        or relative humidity.)
///         Outputs: Rows and columns of the 2x2 matrix and the
///                  reciprocal of the determinant (to calculate inverse):
///                  muRow1Column1
///                  muRow1Column2
///                  muRow2Column1
///                  muRow2Column2
///                  muRecipDeterminant
void eval_moisture_control_traj(atlas::FieldSet & fields);

// TODO(Mayeul) Remove duplicate evalMuThetavTL in control2analysis_linearvarchange.h/cc
/// \details The linear moisture control variable transform
///          that transforms potential temperature increments and
///          total water (qt) increments to moisture control variable
///          and virtual potential temperature increments.
void eval_moisture_control_tl(atlas::FieldSet & incFlds,
                              const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate evalMuThetavAD in control2analysis_linearvarchange.h/cc
/// \details Adjoint of the linear moisture control variable transform
///          that transforms potential temperature increments and
///          total water (qt) increments to moisture control variable
///          and virtual potential temperature increments.
void eval_moisture_control_ad(atlas::FieldSet & hatFlds,
                              const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate evalQtThetaTL in control2analysis_linearvarchange.h/cc
/// \details Inverse of the linear moisture control variable transform.
///          It transforms virtual potential temperature and
///          moisture control variable increments to total water (qt) and
///          potential temperature increments.
void eval_moisture_control_inv_tl(atlas::FieldSet & incFlds,
                                  const atlas::FieldSet & augStateFlds);

// TODO(Mayeul) Remove duplicate evalQtThetaAD in control2analysis_linearvarchange.h/cc
/// \details Adjoint of the inverse of the linear moisture control variable transform
///          that transforms virtual potential temperature and
///          moisture control variable increments to total water (qt) and
///          potential temperature increments.
void eval_moisture_control_inv_ad(atlas::FieldSet & hatFlds,
                                  const atlas::FieldSet & augStateFlds);
}  // namespace mo
