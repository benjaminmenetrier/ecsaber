/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <string>

#include "src/ErrorCovariance.h"
#include "src/Geometry.h"
#include "src/Gom.h"
#include "src/HorizScaleDecomposition.h"
#include "src/Increment.h"
#include "src/IncrEnsCtlVec.h"
#include "src/IncrModCtlVec.h"
#include "src/Interpolator.h"
#include "src/LocalizationMatrix.h"
#include "src/Locations.h"
#include "src/Model.h"
#include "src/ModelBias.h"
#include "src/ModelBiasCovariance.h"
#include "src/ModelBiasCtlVec.h"
#include "src/ModelBiasEstimator.h"
#include "src/ModelBiasIncrement.h"
#include "src/ObsSpace.h"
#include "src/ObsVec.h"
#include "src/State.h"

namespace quench {

struct Traits {
  static std::string name() {return "quench";}
  static std::string nameCovar() {return "quenchCovariance";}

  using Geometry = quench::Geometry;

  using State = quench::State;
  using Model = quench::Model;
  using Increment = quench::Increment;
  using IncrEnsCtlVec = quench::IncrEnsCtlVec;
  using IncrModCtlVec = quench::IncrModCtlVec;
  using Covariance = quench::ErrorCovariance;

  using ModelAuxControl = quench::ModelBias;
  using ModelAuxControlEstimator = quench::ModelBiasEstimator;
  using ModelAuxIncrement = quench::ModelBiasIncrement;
  using ModelAuxCtlVec = quench::ModelBiasCtlVec;
  using ModelAuxCovariance = quench::ModelBiasCovariance;

  using ObsSpace = quench::ObsSpace;
  using ObsVector = quench::ObsVec;

  using GeoVaLs = quench::Gom;
  using Locations = quench::Locations;

  using LocalizationMatrix = quench::LocalizationMatrix;
  using Interpolator = quench::Interpolator;
  using HorizScaleDecomposition = quench::HorizScaleDecomposition;
};

}  // namespace quench