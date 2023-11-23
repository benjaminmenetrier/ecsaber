#!/usr/bin/env bash

# Project name
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Repos
repos="oops
saber
vader"

for repo in ${repos}; do
  # Get hash
  cd ${scriptDir}/../${repo}-jedi
  git rev-parse HEAD > ${scriptDir}/${repo}_hash
done
