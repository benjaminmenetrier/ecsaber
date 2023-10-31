/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/State.h"

#include <vector>

#include "eckit/exception/Exceptions.h"

#include "util/Logger.h"

#include "src/Fields.h"
#include "src/Increment.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const eckit::Configuration & file)
  : fields_()
{
  const oops::Variables vars(file.getStringVector("variables"));
  fields_.reset(new Fields(resol, vars, util::DateTime()));
  if (file.has("filepath")) {
    oops::Log::info() << "Info     : Create state from file" << std::endl;
    fields_->read(file);
  } else {
    oops::Log::info() << "Info     : Create empty state" << std::endl;
    if (file.has("constant value")) {
      fields_->constantValue(file.getDouble("constant value"));
    } else if (file.has("constant group-specific value")) {
      fields_->constantValue(file);
    } else {
      fields_->zero();
    }
  }
  const util::DateTime vt(file.getString("date"));
  fields_->time() = vt;
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const Model &,
             const eckit::Configuration & conf)
  : State(resol, conf) {}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const Tlm &,
             const eckit::Configuration & conf)
  : State(resol, conf) {}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const State & other)
  : fields_(new Fields(*other.fields_, resol))
{
  ASSERT(fields_);
  oops::Log::trace() << "State::State created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const Model &,
             const State & other)
  : State(resol, other) {}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const Model & model,
             const State & other, const eckit::Configuration & conf)
  : State(resol, other) {}
// -----------------------------------------------------------------------------
State::State(const State & other)
  : fields_(new Fields(*other.fields_))
{
  oops::Log::trace() << "State::State copied." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
State & State::operator=(const State & rhs) {
  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
State & State::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  *fields_+=dx.fields();
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void State::read(const eckit::Configuration & files) {
  fields_->read(files);
  fields_->fields().haloExchange();
}
// -----------------------------------------------------------------------------
void State::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
// size_t State::serialSize() const {
//   size_t nn = fields_->serialSize();
//   return nn;
// }
// -----------------------------------------------------------------------------
// void State::serialize(std::vector<double> & vect) const {
//   fields_->serialize(vect);
// }
// -----------------------------------------------------------------------------
// void State::deserialize(const std::vector<double> & vect, size_t & index) {
//   fields_->deserialize(vect, index);
//}
// -----------------------------------------------------------------------------
void State::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << this->validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet accessor
// -----------------------------------------------------------------------------
void State::toFieldSet(atlas::FieldSet & fset) const {
  fields_->toFieldSet(fset);
}
// -----------------------------------------------------------------------------
void State::fromFieldSet(const atlas::FieldSet & fset) {
  fields_->fromFieldSet(fset);
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void State::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void State::accumul(const double & zz, const State & xx) {
  fields_->axpy(zz, *xx.fields_);
}
// -----------------------------------------------------------------------------

}  // namespace quench
