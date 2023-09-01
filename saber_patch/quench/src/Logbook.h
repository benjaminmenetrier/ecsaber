/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/memory/NonCopyable.h"

#include "util/LogbookHelper.h"

namespace quench {

/// Connects to the OOPS execution Logbook
/*!
 * This class allows to connect to oops::LogbookHelper - an execution tracking mechanism;
 * when the execution Logbook is updated, a callback to Logbook update method is invoked,
 * where the Logbook information can be used to update local information; here stored in
 * eckit::LocalConfiguration, but possibly in Fortran for other models (for fast access).
 */
// -----------------------------------------------------------------------------

class Logbook : public util::Printable,
                  private eckit::NonCopyable {
 public:
  static const std::string classname() {return "quench::Logbook";}

/// Initialization, finalization
  static void start();
  static void stop();

/// Destructor
  ~Logbook() {}

 private:
  static Logbook & getInstance();
  Logbook();

/// Basic operators
  void update(const eckit::Configuration & conf);
  void connectToLogbook();

/// Print
  void print(std::ostream & os) const;

/// Data
  std::unique_ptr<eckit::LocalConfiguration> conf_;
};

// -----------------------------------------------------------------------------

Logbook & Logbook::getInstance() {
  static Logbook theLogbook;
  return theLogbook;
}

// -----------------------------------------------------------------------------

Logbook::Logbook() : conf_(new eckit::LocalConfiguration())
{}

// -----------------------------------------------------------------------------

void Logbook::start() {
  oops::Log::info() << "Logbook::strat Logbook starting" << std::endl;
  getInstance().connectToLogbook();
}

// -----------------------------------------------------------------------------

void Logbook::stop() {
  oops::Log::info() << "Logbook::stop Logbook stopping" << std::endl;
}

// -----------------------------------------------------------------------------

void Logbook::update(const eckit::Configuration & conf) {
  conf_.reset(new eckit::LocalConfiguration(conf));
  // MC inspect eckit changes
  // *conf_ = conf;
  oops::Log::info() << "Logbook::update done: " << *this << std::endl;
}

// -----------------------------------------------------------------------------

void Logbook::connectToLogbook() {
  util::LogbookHelper::connectCallback([this]
                                       (const eckit::Configuration & conf)
                                       { return update(conf); });
}

// -----------------------------------------------------------------------------

void Logbook::print(std::ostream & os) const {
  os << *conf_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
