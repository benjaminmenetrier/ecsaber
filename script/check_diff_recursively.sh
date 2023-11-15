#!/usr/bin/env bash

# Input path
projectName=$1
input=$2
commandPath=$3

# Same path for JEDI repo
inputJedi=${input/ecsaber\/${projectName}/ecsaber\/${projectName}-jedi}

# Script directory
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Loop recursively over items
for itemJedi in `find ${inputJedi}`; do
  # Check if item is a file
  if test -f "${itemJedi}"; then
    # Get directory and file
    dirJedi=$(dirname "${itemJedi}")
    file=$(basename "${itemJedi}")

    # Same directory for non-JEDI repo
    dir=${dirJedi/${projectName}-jedi/${projectName}}

    # Update file
    ${scriptDir}/check_diff.sh ${itemJedi} ${dir}/${file} ${commandPath}
  fi
done
