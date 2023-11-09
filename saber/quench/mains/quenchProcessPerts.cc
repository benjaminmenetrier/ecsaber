/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "saber/oops/ProcessPerts.h"
//  #include "src/Logbook.h"  // TODO(Benjamin): only for latest OOPS versions
#include "src/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::ProcessPerts<quench::Traits> processperts;
//  quench::Logbook::start();  // TODO(Benjamin): only for latest OOPS versions
  run.execute(processperts);
//  quench::Logbook::stop();  // TODO(Benjamin): only for latest OOPS versions
  return 0;
}
