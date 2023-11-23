/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <algorithm>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "mo/constants.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/functions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

void eval_sat_vapour_pressure_nl(const std::string & svp_file, atlas::FieldSet & fields)
{
  oops::Log::trace() << "[svp()] starting ..." << std::endl;

  // normalised T lambda function
  // which enforces upper and lower bounds on T
  auto normalisedT = [](const double tVal) {
    double t = (tVal - constants::TLoBound)/constants::Tinc;
    double t1 = (t < 0.0 ? 0.0 : t);
    double t2 = (t1 >= static_cast<double>(constants::svpLookUpLength - 1)) ?
                       static_cast<double>(constants::svpLookUpLength - 1) : t1;
    return t2;
  };

  // lambda function returning the index of a given temperature value;
  // avoids returning last index in table
  auto index = [](const double normalisedTVal) {
    std::size_t i = static_cast<std::size_t>(normalisedTVal);
    return (i == constants::svpLookUpLength - 1 ? constants::svpLookUpLength - 2 : i);
  };

  // weight lamba function
  auto weight = [index](const double normalisedTVal) {
    std::size_t i = index(normalisedTVal);
    return ( normalisedTVal -  static_cast<double>(i) );
  };

  // lambda function to interpolate values from table to given temperature
  auto interp = [](const double weight, const double v1, const double v2) {
     return (weight * v2 + (1 - weight) * v1);
  };

  // lambda function to determine the gradient of variable in table
  // with respect to temperature
  auto normalisedDiff = [](const double v1, const double v2) {
     return (v2 - v1)/ constants::Tinc;
  };

  double w;  // weight
  std::size_t indx;  // index

  const auto tView  = make_view<const double, 2>(fields["air_temperature"]);
  const std::vector<std::string> vars{"svp", "dlsvp", "svpW", "dlsvpW"};
  oops::patch::Variables lookUpVars(vars);
  auto lookUpData = functions::getLookUps(svp_file, lookUpVars,
                               constants::svpLookUpLength);

  const std::vector<std::string> fnames {"svp", "dlsvpdT"};
  int ival = 0;  // set ival = 2 to get svp wrt water
  for (auto & ef : fnames) {
    if (fields.has(ef)) {
      auto svpView = make_view<double, 2>(fields[ef]);

      auto conf = atlas::util::Config("levels", fields[ef].levels()) |
                  atlas::util::Config("include_halo", true);

      // check this recipe to calculate svp is correct
      const std::vector<double> Lookup = lookUpData[ival];
      ++ival;
      auto evaluateSVP = [&] (atlas::idx_t i, atlas::idx_t j) {
      indx = index(normalisedT(tView(i, j)));
      w = weight(normalisedT(tView(i, j)));
      svpView(i, j) = interp(w, Lookup.at(indx), Lookup.at(indx+1)); };

      auto fspace = fields[ef].functionspace();

      functions::parallelFor(fspace, evaluateSVP, conf);
    }
  }

  oops::Log::trace() << "[svp()] ... exit" << std::endl;
}

}  // namespace mo

