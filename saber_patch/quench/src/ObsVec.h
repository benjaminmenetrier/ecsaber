/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {
  class ObsSpace;

// -----------------------------------------------------------------------------
/// ObsVec empty class

class ObsVec : public util::Printable,
               private util::ObjectCounter<ObsVec> {
 public:
  static const std::string classname() {return "quench::ObsVec";}

  explicit ObsVec(const ObsSpace &) {}
  ObsVec(const ObsVec &, const bool copy = true) {}
  ~ObsVec() {}

  ObsVec & operator = (const ObsVec &) {return *this;}
  ObsVec & operator*= (const double &) {return *this;}
  ObsVec & operator+= (const ObsVec &) {return *this;}
  ObsVec & operator-= (const ObsVec &) {return *this;}
  ObsVec & operator*= (const ObsVec &) {return *this;}
  ObsVec & operator/= (const ObsVec &) {return *this;}

  void zero() {}
  void axpy(const double &, const ObsVec &) {}
  void invert() {}
  void random() {}
  double dot_product_with(const ObsVec &) const {return 0.0;}
  double rms() const {return 0.0;}

  unsigned int size() const {return 0;}

// I/O
  void read(const std::string &) {}
  void save(const std::string &) const {}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
