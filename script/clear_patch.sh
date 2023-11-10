#!/usr/bin/env bash

# Script directory
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Go to ECSABER root directory
cd ${scriptDir}/..

for item in `find . -name *.patch`; do
  echo "Remove patch: "${item}
  rm -f ${item}
done
git status
