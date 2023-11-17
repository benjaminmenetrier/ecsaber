/*
 * (C) Crown Copyright 2022 Met Office
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

#include "mo/common_varchange.h"
#include "mo/constants.h"
#include "mo/functions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

bool evalSatVaporPressure(atlas::FieldSet & fields)
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

  // The check for the presence of required input fields will be performed by the Vader
  // algorithm when this code is in a Vader Recipe. At that time this check can be removed.
  if ( !fields.has("air_temperature") ||
       !fields.has("svp")) {
    return false;
  }

  if (fields["svp"].shape(0) != fields["air_temperature"].shape(0)
      || fields["svp"].levels() != fields["air_temperature"].levels()) {
    // svp field not compatible with air temperature field, cannot continue.
    return false;
  }
  const auto tView  = make_view<const double, 2>(fields["air_temperature"]);
  const std::vector<std::string> vars{"svp", "dlsvp", "svpW", "dlsvpW"};
  oops::patch::Variables lookUpVars(vars);
  auto lookUpData = functions::getLookUps(constants::commonVarChangeFilePath, lookUpVars,
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

  return true;
}

bool evalSatSpecificHumidity(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[getQsat()] starting ..." << std::endl;

  const auto pbarView = make_view<const double, 2>(fields["air_pressure"]);
  const auto svpView = make_view<const double, 2>(fields["svp"]);
  const auto tView = make_view<const double, 2>(fields["air_temperature"]);
  auto qsatView = make_view<double, 2>(fields["qsat"]);

  auto conf = atlas::util::Config("levels", fields["qsat"].levels()) |
              atlas::util::Config("include_halo", true);

  double fsubw;
  double svp_times_fsubw;
  auto evaluateQsat = [&] (atlas::idx_t i, atlas::idx_t j) {
    // Compute the factor that converts from sat vapour pressure in a
    // pure water system to sat vapour pressure in air, FSUBW. This formula
    // is taken from equation A4.7 of Adrian Gill's book: Atmosphere-Ocean
    // Dynamics. Note that his formula works in terms of pressure in mb and
    // temperature in Celsius, so conversion of units leads to the slightly
    // different equation used here.
    fsubw = 1.0 + 1.0E-8 * pbarView(i, j) * (4.5 +
            6.0e-4 * (tView(i, j) - constants::zerodegc) *
                     (tView(i, j) - constants::zerodegc));

    // Multiply by FSUBW to convert to saturated vapour pressure in air
    // (equation A4.6 of Adrian Gill's book)
    svp_times_fsubw = fsubw * svpView(i, j);

    // Now form the accurate expression for QS, which is a rearranged
    // version of equation A4.3 of Gill's book.
    // Note that at very low pressures we apply a fix, to prevent a
    // singularity (Qsat tends to 1.0 kg/kg).
    qsatView(i, j) = constants::rd_over_rv * svp_times_fsubw /
          (std::max(pbarView(i, j), svp_times_fsubw) -
          (1.0 - constants::rd_over_rv) * svp_times_fsubw);
  };

  auto fspace = fields["qsat"].functionspace();

  functions::parallelFor(fspace, evaluateQsat, conf);

  oops::Log::trace() << "[getQsat()] ... exit" << std::endl;

  return true;
}

bool evalAirPressureLevels(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalAirPressureLevels()] starting ..." << std::endl;

  const auto ds_elmo = make_view<const double, 2>(fields["exner_levels_minus_one"]);
  const auto ds_plmo = make_view<const double, 2>(fields["air_pressure_levels_minus_one"]);
  const auto ds_t = make_view<const double, 2>(fields["potential_temperature"]);
  const auto ds_hl = make_view<const double, 2>(fields["height_levels"]);
  auto ds_pl = make_view<double, 2>(fields["air_pressure_levels"]);

  idx_t levels(fields["air_pressure_levels"].levels());
  for (idx_t jn = 0; jn < fields["air_pressure_levels"].shape(0); ++jn) {
    for (idx_t jl = 0; jl < levels - 1; ++jl) {
      ds_pl(jn, jl) = ds_plmo(jn, jl);
    }

    // Note that I am calculating the exner pressure above the top first and then
    // converting it to pressure
    // Note that strictly speaking we should be using virtual potential temperature here
    // but given that there should be no moisture at the top of the model, we should be
    // able to ignore that contribution.
    //
    // exner^k+1 = exner^k - g (height^k+1 - height^k) / (cp theta_v)
    //
    // pressure^k+1 = reference_pressure * (exner^k+1)**((1.0 / constants::rd_over_cp)
    //
    // where k is the model level index on half levels just below model top.

    ds_pl(jn, levels-1) =  constants::p_zero * pow(
      ds_elmo(jn, levels-2) - (constants::grav * (ds_hl(jn, levels-1) - ds_hl(jn, levels-2))) /
      (constants::cp * ds_t(jn, levels-2)), (1.0 / constants::rd_over_cp));

    ds_pl(jn, levels-1) = ds_pl(jn, levels-1) > 0.0 ? ds_pl(jn, levels-1) : constants::deps;
  }

  oops::Log::trace() << "[evalAirPressureLevels()] ... exit" << std::endl;

  return true;
}

}  // namespace mo
