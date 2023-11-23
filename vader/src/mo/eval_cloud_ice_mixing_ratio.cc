/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "mo/eval_cloud_ice_mixing_ratio.h"
#include "mo/eval_ratio.h"

#include "oops/util/Logger.h"

namespace mo {

// --------------------------------------------------------------------------------------

bool eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
        atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl()] starting ..."
          << std::endl;

  std::vector<std::string> fnames {"m_ci", "m_t",
                                   "mass_content_of_cloud_ice_in_atmosphere_layer"};

  const bool rvalue = eval_ratio_to_second(stateFlds, fnames);

  oops::Log::trace()
          << "[eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl()] ... exit"
          << std::endl;
  return rvalue;
}

// --------------------------------------------------------------------------------------

void eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_tl(
        atlas::FieldSet & incFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_tl()] starting ..."
          << std::endl;

  std::vector<std::string> fnames {"m_ci", "m_t",
                                   "mass_content_of_cloud_ice_in_atmosphere_layer"};

  eval_ratio_to_second_tl(incFlds, stateFlds, fnames);

  oops::Log::trace()
          << "[eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_tl()] ... exit"
          << std::endl;
}
// --------------------------------------------------------------------------------------

void eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_ad(
        atlas::FieldSet & hatFlds,
        const atlas::FieldSet & stateFlds) {
  oops::Log::trace()
          << "[eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_ad()] starting ..."
          << std::endl;

  std::vector<std::string> fnames {"m_ci", "m_t",
                                   "mass_content_of_cloud_ice_in_atmosphere_layer"};

  eval_ratio_to_second_ad(hatFlds, stateFlds, fnames);

  oops::Log::trace()
          << "[eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_ad()] ... exit"
          << std::endl;
}

// --------------------------------------------------------------------------------------

}  // namespace mo
