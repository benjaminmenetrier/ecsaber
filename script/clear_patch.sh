#!/usr/bin/env bash

# Script directory
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Go to ECSABER root directory
cd ${scriptDir}/..

for item in `find . -name *.patch`; do
  echo "Remove patch file: "${item}
  rm -f ${item}
done

for item in `find . -name *.tmp`; do
  echo "Remove temporary file: "${tmp}
  rm -f ${item}
done

for item in `find . -name *.bak`; do
  echo "Remove backup file: "${item}
  rm -f ${item}
done
