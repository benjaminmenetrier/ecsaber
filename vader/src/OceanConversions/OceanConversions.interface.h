/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_OCEANCONVERSIONS_OCEANCONVERSIONS_INTERFACE_H_
#define SRC_OCEANCONVERSIONS_OCEANCONVERSIONS_INTERFACE_H_

namespace vader {

extern "C" {
// call Fortran code
  const float gsw_rho_t_exact_f90(const float & sal,
                                  const float & temp,
                                  const float & pressure);
  const float gsw_p_from_z_f90(const float & depth,
                               const float & latitude);
  const float gsw_pt_from_t_f90(const float & sal,
                                const float & temp,
                                const float & pressure);
  const float gsw_ct_from_t_f90(const float & sal,
                                const float & temp,
                                const float & pressure);
  const float gsw_ct_from_pt_f90(const float & sal,
                                const float & temp);
  const float gsw_t_from_ct_f90(const float & sal,
                                const float & temp,
                                const float & pressure);
  const float gsw_sa_from_sp_f90(const float & sal,
                                 const float & pressure,
                                 const float & longitude,
                                 const float & latitude);
}  // extern C

}  // namespace vader

#endif  // SRC_OCEANCONVERSIONS_OCEANCONVERSIONS_INTERFACE_H_
