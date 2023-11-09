/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "src/Fields.h"

#include "util/abor1_cpp.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {
  class Gom;
  class Geometry;
  class Increment;
  class Locations;
  class Model;
  class Tlm;
  class Variables;

///  model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class State : public util::Printable,
              private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "quench::State";}

/// Constructor, destructor
  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const Model &,
          const eckit::Configuration &);
  State(const Geometry &, const Tlm &,
          const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const Geometry &, const Model &,
          const State &);
  State(const Geometry &, const Model &,
          const State &, const eckit::Configuration &);
  State(const State &);
  virtual ~State() {}
  State & operator=(const State &);

/// Interpolate to observation location
  void interpolate(const Locations &, Gom &) const
    {ABORT("not implemented yet");}
/// Interpolate full fields
  void forceWith(const State &, const Variables &)
    {ABORT("not implemented yet");}

/// Interactions with Increment
  State & operator+=(const Increment &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const {
    return fields_->geometry();
  }

/// Other
  void activateModel()
    {ABORT("not implemented yet");}
  void deactivateModel()
    {ABORT("not implemented yet");}

  void zero();
  void accumul(const double &, const State &);

/// ATLAS FieldSet accessor
  const atlas::FieldSet & fieldSet() const {return fields_->fieldSet();}
  atlas::FieldSet & fieldSet() {return fields_->fieldSet();}
  void synchronizeFields() {fields_->synchronizeFields();}

 private:
  std::unique_ptr<Fields> fields_;
  void print(std::ostream &) const;
};
// -----------------------------------------------------------------------------

}  // namespace quench
