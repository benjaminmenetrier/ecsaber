/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/serialisation/Stream.h"

#include "src/Variables.h"

#include "util/abor1_cpp.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class Geometry;

// -----------------------------------------------------------------------------
/// Class to represent a Fields for the  model
class Fields : public util::Printable,
               private util::ObjectCounter<Fields> {
 public:
  static const std::string classname() {return "quench::Fields";}

// Constructors
  Fields(const Geometry &, const Variables &, const util::DateTime &);
  Fields(const Fields &, const Geometry &);
  Fields(const Fields &, const bool);
  Fields(const Fields &);
  ~Fields() {}

// Basic operators
  void zero();
  void constantValue(const double &);
  Fields & operator=(const Fields &);
  Fields & operator+=(const Fields &);
  Fields & operator-=(const Fields &);
  Fields & operator*=(const double &);
  void axpy(const double &, const Fields &);
  double dot_product_with(const Fields &) const;
  void schur_product_with(const Fields &);
  void dirac(const eckit::Configuration &);
  void random();
  void diff(const Fields &, const Fields &);

  friend eckit::Stream & operator<<(eckit::Stream &, const Fields &);
  friend eckit::Stream & operator>>(eckit::Stream &, Fields &);

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  double min(const Variables & var) const;
  double max(const Variables & var) const;
  std::shared_ptr<const Geometry> geometry() const {return geom_;}
  const Variables & variables() const {return vars_;}
  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}

/// Serialization
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

/// ATLAS FieldSet
  const atlas::FieldSet & fieldSet() const {return fset_;}
  atlas::FieldSet & fieldSet() {return fset_;}
  void synchronizeFields();

 private:
  void print(std::ostream &) const;
  std::shared_ptr<const Geometry> geom_;
  const Variables vars_;
  util::DateTime time_;
  mutable atlas::FieldSet fset_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
