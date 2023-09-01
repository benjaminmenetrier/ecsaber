#!/usr/bin/env bash
# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Create logbook
if [ -f logbook.json ]; then
  rm logbook.json
fi

# Set number of OpenMP threads
export OMP_NUM_THREADS=$3

# Log and test outputs extensions
flog=$5.log.out
ftest=$5.test.out

if test -f "$5"; then
  # Run job, create log and test outputs
  mpirun -n $2 $4 | tee ${flog} && \
  grep 'Test     : ' ${flog} > ${ftest}

  # Compare reference and log output
  $1 ${flog} $5 1.0e-12 0
else
  # Run job
  mpirun -n $2 $4 | tee ${flog}
fi
