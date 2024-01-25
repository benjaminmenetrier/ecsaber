
/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace.h"

#include "mo/constants.h"
#include "mo/functions.h"
#include "mo/model2geovals_varchange.h"

#include "oops/util/Logger.h"

using atlas::array::make_view;
using atlas::idx_t;
using atlas::util::Config;

namespace mo {

void initField_rank2(atlas::Field & field, const double value_init)
{
  setUniformValue_rank2(field, value_init);
}


void setUniformValue_rank2(atlas::Field & field, const double value)
{
  auto ds_view = make_view<double, 2>(field);

  ds_view.assign(value);
}


bool evalTotalMassMoistAir(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalTotalMassMoistAir()] starting ..." << std::endl;

  const auto ds_m_v  = make_view<const double, 2>(fields["m_v"]);
  const auto ds_m_ci = make_view<const double, 2>(fields["m_ci"]);
  const auto ds_m_cl = make_view<const double, 2>(fields["m_cl"]);
  const auto ds_m_r  = make_view<const double, 2>(fields["m_r"]);
  auto ds_m_t  = make_view<double, 2>(fields["m_t"]);

  auto fspace = fields["m_t"].functionspace();

  auto evaluateMt = [&] (idx_t i, idx_t j) {
    ds_m_t(i, j) = 1 + ds_m_v(i, j) + ds_m_ci(i, j) + ds_m_cl(i, j) + ds_m_r(i, j); };

  auto conf = Config("levels", fields["m_t"].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateMt, conf);

  oops::Log::trace() << "[evalTotalMassMoistAir()] ... exit" << std::endl;

  return true;
}

bool evalRatioToMt(atlas::FieldSet & fields, const std::vector<std::string> & vars)
{
  oops::Log::trace() << "[evalRatioToMt()] starting ..." << std::endl;

  // fields[0] = m_x = [ mv | mci | mcl | m_r ]
  const auto ds_m_x  = make_view<const double, 2>(fields[vars[0]]);
  const auto ds_m_t  = make_view<const double, 2>(fields[vars[1]]);
  auto ds_tfield  = make_view<double, 2>(fields[vars[2]]);

  auto fspace = fields[vars[2]].functionspace();

  auto evaluateRatioToMt = [&] (idx_t i, idx_t j) {
    ds_tfield(i, j) = ds_m_x(i, j) / ds_m_t(i, j);
  };

  auto conf = Config("levels", fields[vars[2]].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateRatioToMt, conf);

  oops::Log::trace() << "[evalRatioToMt()] ... exit" << std::endl;

  return true;
}


bool evalSpecificHumidity(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalSpecificHumidity()] starting ..." << std::endl;

  std::vector<std::string> fnames {"m_v", "m_t", "specific_humidity"};

  bool rvalue = evalRatioToMt(fields, fnames);

  oops::Log::trace() << "[evalSpecificHumidity()] ... exit" << std::endl;

  return rvalue;
}


bool evalTotalRelativeHumidity(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalTotalRelativeHumidity()] starting ..." << std::endl;

  const auto qView = make_view<const double, 2>(fields["specific_humidity"]);
  const auto qclView = make_view<const double, 2>
                 (fields["mass_content_of_cloud_liquid_water_in_atmosphere_layer"]);
  const auto qciView = make_view<const double, 2>
                 (fields["mass_content_of_cloud_ice_in_atmosphere_layer"]);
  const auto qrainView = make_view<const double, 2>(fields["qrain"]);
  const auto qsatView = make_view<const double, 2>(fields["qsat"]);
  auto rhtView = make_view<double, 2>(fields["rht"]);

  auto conf = Config("levels", fields["rht"].levels()) |
              Config("include_halo", true);

  auto evaluateRHT = [&] (idx_t i, idx_t j) {
    rhtView(i, j) = (qView(i, j) + qclView(i, j) + qciView(i, j)
    + qrainView(i, j)) / qsatView(i, j) * 100.0;

    if (rhtView(i, j) < 0.0)
    {
        rhtView(i, j) = 0.0;
    }
  };

  auto fspace = fields["rht"].functionspace();

  functions::parallelFor(fspace, evaluateRHT, conf);

  oops::Log::trace() << "[evalTotalRelativeHumidity()] ... exit" << std::endl;

  return true;
}

bool evalMassCloudIce(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalMassCloudIce()] starting ..." << std::endl;

  std::vector<std::string> fnames {"m_ci", "m_t",
                                   "mass_content_of_cloud_ice_in_atmosphere_layer"};
  bool rvalue = evalRatioToMt(fields, fnames);

  oops::Log::trace() << "[evalMassCloudIce()] ... exit" << std::endl;

  return rvalue;
}


bool evalMassCloudLiquid(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalMassCloudLiquid()] starting ..." << std::endl;

  std::vector<std::string> fnames {"m_cl", "m_t",
                                   "mass_content_of_cloud_liquid_water_in_atmosphere_layer"};
  bool rvalue = evalRatioToMt(fields, fnames);

  oops::Log::trace() << "[evalMassCloudLiquid()] ... exit" << std::endl;

  return rvalue;
}


bool evalMassRain(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalMassRain()] starting ..." << std::endl;

  std::vector<std::string> fnames {"m_r", "m_t", "qrain"};
  bool rvalue = evalRatioToMt(fields, fnames);

  oops::Log::trace() << "[evalMassRain()] ... exit" << std::endl;

  return rvalue;
}


bool evalAirTemperature(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalAirTemperature()] starting ..." << std::endl;

  const auto ds_theta  = make_view<const double, 2>(fields["potential_temperature"]);
  const auto ds_exner = make_view<const double, 2>(fields["exner"]);
  auto ds_atemp = make_view<double, 2>(fields["air_temperature"]);

  auto fspace = fields["air_temperature"].functionspace();

  auto evaluateAirTemp = [&] (idx_t i, idx_t j) {
    ds_atemp(i, j) = ds_theta(i, j) * ds_exner(i, j); };

  auto conf = Config("levels", fields["air_temperature"].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateAirTemp, conf);

  oops::Log::trace() << "[evalAirTemperature()] ... exit" << std::endl;

  return true;
}

bool evalSpecificHumidityFromRH_2m(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalSpecificHumidityFromRH_2m()] starting ..." << std::endl;

  const auto ds_qsat = make_view<const double, 2>(fields["qsat"]);
  const auto ds_rh = make_view<const double, 2>(fields["relative_humidity_2m"]);
  auto ds_q2m = make_view<double, 2>(fields["specific_humidity_at_two_meters_above_surface"]);

  auto fspace = fields["specific_humidity_at_two_meters_above_surface"].functionspace();

  auto evaluateSpecificHumidity_2m = [&] (idx_t i, idx_t j) {
    ds_q2m(i, j) = 0.01 * ds_rh(i, j) * ds_qsat(i, j); };

  auto conf = Config("levels",
    fields["specific_humidity_at_two_meters_above_surface"].levels()) |
              Config("include_halo", true);

  functions::parallelFor(fspace, evaluateSpecificHumidity_2m, conf);

  oops::Log::trace() << "[evalSpecificHumidityFromRH_2m()] ... exit" << std::endl;

  return true;
}


bool evalParamAParamB(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalParamAParamB2()] starting ..." << std::endl;

  std::size_t blindex;
  if (!fields["height"].metadata().has("boundary_layer_index")) {
    oops::Log::error() << "ERROR - data validation failed "
                          "we expect boundary_layer_index value "
                          "in the meta data of the height field" << std::endl;
  }
  fields["height"].metadata().get("boundary_layer_index", blindex);

  const auto heightView = make_view<const double, 2>(fields["height"]);
  const auto heightLevelsView = make_view<const double, 2>(fields["height_levels"]);
  const auto pressureLevelsView = make_view<const double, 2>
      (fields["air_pressure_levels_minus_one"]);
  const auto specificHumidityView = make_view<const double, 2>(fields["specific_humidity"]);
  auto param_aView = make_view<double, 2>(fields["surf_param_a"]);
  auto param_bView = make_view<double, 2>(fields["surf_param_b"]);

  // temperature at level above boundary layer
  double t_bl;
  // temperature at model surface height
  double t_msh;

  double exp_pmsh = constants::Lclr * constants::rd / constants::grav;

  for (idx_t jn = 0; jn < param_aView.shape(0); ++jn) {
    t_bl = (-constants::grav / constants::rd) *
           (heightLevelsView(jn, blindex + 1) - heightLevelsView(jn, blindex)) /
           log(pressureLevelsView(jn, blindex + 1) / pressureLevelsView(jn, blindex));

    t_bl = t_bl / (1.0 + constants::c_virtual * specificHumidityView(jn, blindex));

    t_msh = t_bl + constants::Lclr * (heightView(jn, blindex) - heightLevelsView(jn, 0));

    param_aView(jn, 0) = heightLevelsView(jn, 0) + t_msh / constants::Lclr;
    param_bView(jn, 0) = t_msh / (pow(pressureLevelsView(jn, 0), exp_pmsh) * constants::Lclr);
  }

  oops::Log::trace() << "[evalParamAParamB()] ... exit" << std::endl;

  return true;
}

bool evalGpHeightMinusOne(atlas::FieldSet & fields)
{
  oops::Log::trace() << "[evalGpHeightMinusOne()] starting ..." << std::endl;

  const auto height_lmoView = make_view<const double, 2>(fields["height_levels_minus_one"]);
  auto gph_lmoView = make_view<double, 2>(fields["geopotential_height_levels_minus_one"]);

  gph_lmoView.assign(height_lmoView);

  oops::Log::trace() << "[evalGpHeightMinusOne()] ... exit" << std::endl;

  return true;
}

}  // namespace mo
