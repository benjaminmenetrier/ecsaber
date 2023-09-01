/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "oops/runs/ConvertState.h"
#include "oops/runs/Run.h"
//  #include "src/Logbook.h"  // TODO(Benjamin): only for latest OOPS versions
#include "src/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::ConvertState<quench::Traits> cs;
//  quench::Logbook::start();  // TODO(Benjamin): only for latest OOPS versions
  run.execute(cs);
//  quench::Logbook::stop();  // TODO(Benjamin): only for latest OOPS versions
  return 0;
}
