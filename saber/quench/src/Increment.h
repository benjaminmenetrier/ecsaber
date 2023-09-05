/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>

#include "eckit/serialisation/Stream.h"

#include "oops/base/GeneralizedDepartures.h"

#include "util/abor1_cpp.h"
#include "util/DateTime.h"
#include "util/dot_product.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "src/Fields.h"
#include "src/Geometry.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {
  class Gom;
  class Locations;
  class ModelBiasIncrement;
  class ErrorCovariance;
  class State;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class Increment : public oops::GeneralizedDepartures,
                  public util::Printable,
                  private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "quench::Increment";}

/// Constructor, destructor
  Increment(const Geometry &, const oops::Variables &, const util::DateTime &);
  Increment(const Geometry &, const oops::Variables &, const util::DateTime &,
            const util::DateTime &);
  Increment(const Geometry &, const Increment &);
  Increment(const Increment &, const bool);
  Increment(const Increment &) {}
  virtual ~Increment() {}

/// Basic operators
  void diff(const State &, const State &);
  void zero();
  void zero(const util::DateTime &);
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true);
  double dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();

/// Interpolate to observation location
  void interpolateTL(const Locations &, Gom &) const
    {ABORT("not implemented yet");}
  void interpolateAD(const Locations &, const Gom &)
    {ABORT("not implemented yet");}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  double max(const oops::Variables & var) const {return fields_->max(var);}
  double min(const oops::Variables & var) const {return fields_->min(var);}
  void dirac(const eckit::Configuration &);
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const {
    return fields_->geometry();
  }

/// Serialization
  friend eckit::Stream & operator<<(eckit::Stream & s, const Increment & dx)
    {s << dx.fields(); return s;}
  friend eckit::Stream & operator>>(eckit::Stream & s, Increment & dx)
    {s >> dx.fields(); return s;}

/// Other
  void activateModel()
    {ABORT("not implemented yet");}
  void deactivateModel()
    {ABORT("not implemented yet");}

  void accumul(const double &, const State &);

/// ATLAS FieldSet accessor
  void toFieldSet(atlas::FieldSet &) const;
  void toFieldSetAD(const atlas::FieldSet &) {ABORT("not implemented");}
  void fromFieldSet(const atlas::FieldSet &);

/// Data
 private:
  std::unique_ptr<Fields> fields_;
  void print(std::ostream &) const;
};
// -----------------------------------------------------------------------------

}  // namespace quench
