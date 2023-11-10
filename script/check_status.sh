#!/usr/bin/env bash

# Project name
projectName=$1

# Script directory
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Reference hash
refHash=`cat ${scriptDir}/${projectName}_hash`

# Go to ${projectName}-jedi
cd ${scriptDir}/../${projectName}-jedi

# Initialize ${projectName}-jedi status
correct=true

# Check for non-committed modifications
modifiedFiles=`git diff --name-only | wc -l`
if test "${modifiedFiles}" != "0"; then
  correct=false
else
  # Check for git hash
  currentHash=`git rev-parse HEAD`
  if test "${currentHash}" != "${refHash}"; then
    correct=false
  fi
fi

## Return status
if test "${correct}" = "true"; then
  echo -e "-- ${projectName}-jedi status is OK, update not required"
  exit 0
else
  echo -e "-- WARNING: ${projectName}-jedi status is not OK, update required"
  exit 1
fi
