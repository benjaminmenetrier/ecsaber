/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/array.h"
#include "atlas/parallel/omp/omp.h"

#include "mo/constants.h"
#include "mo/eval_moisture_control.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace mo {

void eval_moisture_control_traj(atlas::FieldSet & fields) {
  auto qtView = make_view<const double, 2>(fields["qt"]);
  auto qView = make_view<const double, 2>(fields["specific_humidity"]);
  auto thetaView = make_view<const double, 2>(fields["potential_temperature"]);
  auto exnerView = make_view<const double, 2>(fields["exner"]);
  auto dlsvpdTView = make_view<const double, 2>(fields["dlsvpdT"]);
  auto qsatView = make_view<const double, 2>(fields["qsat"]);
  auto muAView = make_view<const double, 2>(fields["muA"]);
  auto muH1View = make_view<const double, 2>(fields["muH1"]);

  // muRow1Column1, muRow1Column2, muRow2Column1, muRow2Column2
  // are coefficients of a (2x2) matrix = A
  //  (mu')       = A (qt')     where A is
  //  (theta_v')      (theta')
  //
  //  ( muA/qsat    - (muA/qsat) muH1 qT exner_bar dlsvpdT )
  //  (                                                    )
  //  (c_v theta     (1 + cv q)                            )
  //
  auto muRow1Column1View = make_view<double, 2>(fields["muRow1Column1"]);
  auto muRow1Column2View = make_view<double, 2>(fields["muRow1Column2"]);
  auto muRow2Column1View = make_view<double, 2>(fields["muRow2Column1"]);
  auto muRow2Column2View = make_view<double, 2>(fields["muRow2Column2"]);
  auto muRecipDeterminantView = make_view<double, 2>(fields["muRecipDeterminant"]);

  // the comments below are there to allow checking with the VAR code.
  const atlas::idx_t n_levels(fields["potential_temperature"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < thetaView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      muRow1Column1View(ih, ilev) = muAView(ih, ilev) / qsatView(ih, ilev);  // beta2 * muA
      muRow1Column2View(ih, ilev) = - qtView(ih, ilev)  * muH1View(ih, ilev)
        * exnerView(ih, ilev) * dlsvpdTView(ih, ilev) * muRow1Column1View(ih, ilev);
      // alpha2 * muA
      muRow2Column1View(ih, ilev) = constants::c_virtual * thetaView(ih, ilev);   // beta1
      muRow2Column2View(ih, ilev) = 1.0 + constants::c_virtual * qView(ih, ilev);  // alpha1
      muRecipDeterminantView(ih, ilev) = 1.0 /(
        muRow2Column2View(ih, ilev) * muRow1Column1View(ih, ilev)
        - muRow1Column2View(ih, ilev) * muRow2Column1View(ih, ilev));
           // 1/( alpha1 * beta2 * muA - alpha2 * muA * beta1)
    }
  }
}

void eval_moisture_control_tl(atlas::FieldSet & incFlds,
                              const atlas::FieldSet & augStateFlds) {
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto thetaIncView = make_view<const double, 2>(incFlds["potential_temperature"]);
  auto qtIncView = make_view<const double, 2>(incFlds["qt"]);
  auto muIncView = make_view<double, 2>(incFlds["mu"]);
  auto thetavIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  const atlas::idx_t n_levels(incFlds["mu"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < muIncView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      muIncView(ih, ilev) = muRow1Column1View(ih, ilev) * qtIncView(ih, ilev)
                        + muRow1Column2View(ih, ilev) * thetaIncView(ih, ilev);
      thetavIncView(ih, ilev) = muRow2Column1View(ih, ilev) * qtIncView(ih, ilev)
                            + muRow2Column2View(ih, ilev) * thetaIncView(ih, ilev);
    }
  }
}

void eval_moisture_control_ad(atlas::FieldSet & hatFlds,
                              const atlas::FieldSet & augStateFlds) {
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);
  auto qtHatView = make_view<double, 2>(hatFlds["qt"]);
  auto muHatView = make_view<double, 2>(hatFlds["mu"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);

  const atlas::idx_t n_levels(hatFlds["mu"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < muHatView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      thetaHatView(ih, ilev) += muRow2Column2View(ih, ilev) * thetavHatView(ih, ilev);
      qtHatView(ih, ilev) += muRow2Column1View(ih, ilev) * thetavHatView(ih, ilev);
      thetaHatView(ih, ilev) += muRow1Column2View(ih, ilev) * muHatView(ih, ilev);
      qtHatView(ih, ilev) += muRow1Column1View(ih, ilev) * muHatView(ih, ilev);
      thetavHatView(ih, ilev) = 0.0;
      muHatView(ih, ilev) = 0.0;
    }
  }
}

void eval_moisture_control_inv_tl(atlas::FieldSet & incFlds,
                                  const atlas::FieldSet & augStateFlds) {
  // Using Cramer's rule to calculate inverse.
  auto muRecipDeterView = make_view<const double, 2>(augStateFlds["muRecipDeterminant"]);
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View  = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto muIncView = make_view<const double, 2>(incFlds["mu"]);
  auto thetavIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto qtIncView = make_view<double, 2>(incFlds["qt"]);
  auto thetaIncView = make_view<double, 2>(incFlds["potential_temperature"]);

  const atlas::idx_t n_levels(incFlds["mu"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < muIncView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      // VAR equivalent in Var_UpPFtheta_qT.f90 for thetaIncView
      // (beta2 * muA * theta_v' +   beta1 * mu') /
      // (alpha1 * beta2 * muA - alpha2 * muA * beta1)
      thetaIncView(ih, ilev) =  muRecipDeterView(ih, ilev) * (
                             muRow1Column1View(ih, ilev) * thetavIncView(ih, ilev)
                           - muRow2Column1View(ih, ilev) * muIncView(ih, ilev) );

      // VAR equivalent in Var_UpPFtheta_qT.f90 for qtIncView
      // (alpha1 * mu_v' -   alpha2 * muA * thetav') /
      // (alpha1 * beta2 * muA - alpha2 * muA * beta1)
      qtIncView(ih, ilev) =  muRecipDeterView(ih, ilev) * (
                           muRow2Column2View(ih, ilev) * muIncView(ih, ilev) -
                           muRow1Column2View(ih, ilev) * thetavIncView(ih, ilev) );
    }
  }
}

void eval_moisture_control_inv_ad(atlas::FieldSet & hatFlds,
                                  const atlas::FieldSet & augStateFlds) {
  auto muRecipDeterView = make_view<const double, 2>(augStateFlds["muRecipDeterminant"]);
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View  = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto qtHatView = make_view<double, 2>(hatFlds["qt"]);
  auto muHatView = make_view<double, 2>(hatFlds["mu"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["potential_temperature"]);

  const atlas::idx_t n_levels(hatFlds["mu"].levels());
  atlas_omp_parallel_for(idx_t ih = 0; ih < muHatView.shape(0); ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      thetavHatView(ih, ilev) += muRecipDeterView(ih, ilev) *
                               muRow1Column1View(ih, ilev) * thetaHatView(ih, ilev);
      muHatView(ih, ilev) -= muRecipDeterView(ih, ilev) *
                           muRow2Column1View(ih, ilev) * thetaHatView(ih, ilev);
      thetavHatView(ih, ilev) -= muRecipDeterView(ih, ilev) *
                               muRow1Column2View(ih, ilev) * qtHatView(ih, ilev);
      muHatView(ih, ilev) += muRecipDeterView(ih, ilev) *
                           muRow2Column2View(ih, ilev) * qtHatView(ih, ilev);
      thetaHatView(ih, ilev) = 0.0;
      qtHatView(ih, ilev) = 0.0;
    }
  }
}

}  // namespace mo
