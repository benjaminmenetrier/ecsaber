/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/assimilation/State4D.h"
#include "oops/base/DataSetBase.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"

namespace oops {

// -----------------------------------------------------------------------------
class FieldSets : public DataSetBase<FieldSet3D, atlas::FunctionSpace> {
  typedef DataSetBase<FieldSet3D, atlas::FunctionSpace> Base_;

 public:
  /// @brief Creates a FieldSets from the Ensemble.
  template<typename MODEL> FieldSets(const State4D<MODEL> &,
                                     const std::vector<int> &);

  FieldSets(const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
                     const std::vector<int> &, const eckit::mpi::Comm &);

  /// @brief  Multiplies each FieldSet3D in this FieldSets with the \p other.
  FieldSets & operator*=(const oops::FieldSet3D & other);
  FieldSets & operator*=(const double zz);

 private:
  std::string classname() const {return "FieldSets";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSets::FieldSets(const State4D<MODEL> & xb, const std::vector<int> & ensmems)
  : Base_(xb.times(), eckit::mpi::self(), ensmems, eckit::mpi::self()) {
  for (size_t jj = 0; jj < ensmems.size(); ++jj) {
    for (unsigned jsub = 0; jsub < xb.times().size(); ++jsub) {
      this->dataset().emplace_back(std::make_unique<FieldSet3D>(xb[jsub].validTime(),
                                                  eckit::mpi::comm()));
    }
  }
  for (unsigned jsub = 0; jsub < xb.times().size(); ++jsub) {
    std::shared_ptr<Ensemble<MODEL>> ensemble =
      EnsemblesCollection<MODEL>::getInstance()[xb[jsub].validTime()];
    ASSERT(ensemble->size() == ensmems.size());
    for (size_t jj = 0; jj < ensemble->size(); ++jj) {
      size_t index = jj*xb.times().size()+jsub;
      this->dataset()[index]->shallowCopy((*ensemble)[jj].increment().fieldSet());
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
