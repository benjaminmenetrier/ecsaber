/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace mo {
namespace functions {
//--
// ++ Atlas Function Spaces ++

/// \brief procedure to call a functor for a given concrete implementation
///        of a function space type
template<typename Functor>
void executeFunc(const atlas::FunctionSpace & fspace, const Functor & functor) {
  if (atlas::functionspace::NodeColumns(fspace)) {
    functor(atlas::functionspace::CubedSphereNodeColumns(fspace));
  } else if (atlas::functionspace::CellColumns(fspace)) {
    functor(atlas::functionspace::CubedSphereCellColumns(fspace));
  } else if (atlas::functionspace::StructuredColumns(fspace)) {
    functor(atlas::functionspace::StructuredColumns(fspace));
  } else {
    oops::Log::error() << "ERROR - a functor call failed "
                          "(function space type not allowed)" << std::endl;
    throw std::runtime_error("a functor call failed");
  }
}

/// \brief wrapper for 'parallel_for'
template<typename Functor>
void parallelFor(const atlas::FunctionSpace & fspace,
                 const Functor& functor,
                 const atlas::util::Config& conf = atlas::util::Config()) {
  executeFunc(fspace, [&](const auto& fspace){fspace.parallel_for(conf, functor);});
}


//--
// ++ I/O processing ++

/// \brief function to read data from a netcdf file
/// sVPFilePath: the path and name of the netcdf file
/// shortname: array to be read from the file
/// lookupSize: dimension of the array
///
std::vector<double> getLookUp(const std::string & sVPFilePath,
                              const std::string & shortName,
                              const std::size_t lookupSize);

/// \brief function to read data from a netcdf file
/// sVPFilePath: the path and name of the netcdf file
/// shortname: list of arrays, with same dimensions, to be read from the file
/// lookupSize: dimension of each array
///
std::vector<std::vector<double>> getLookUps(const std::string & sVPFilePath,
                                            const oops::patch::Variables & vars,
                                            const std::size_t lookupSize);

/// \details getMIOFields returns the effective cloud fractions
///          for the moisture incrementing operator (MIO)
void getMIOFields(atlas::FieldSet & augStateFlds);

/// \details This extracts the scaling coefficients that are applied to Cleff and Cfeff
///          to generate the qcl and qcf increments in the moisture incrementing operator (MIO)
///          The string s can be "qcl_coef" or "qcf_coef"
Eigen::MatrixXd createMIOCoeff(const std::string mioFileName,
                               const std::string s);

extern "C" {
  void umGetLookUp_f90(
    const int &,
    const char *,
    const int &,
    const char *,
    const int &,
    double &);

  void umGetLookUp2D_f90(
    const int &,
    const char *,
    const int &,
    const char *,
    const int &,
    const int &,
    double &);
}  // extern "C"

}  // namespace functions
}  // namespace mo
