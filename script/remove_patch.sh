#!/usr/bin/env bash

# Script directory
scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

find ${inputJedi}



# Source and destination paths
srcPath=$1
dstPath=$2

# Initial check
if test -f "${dstPath}"; then
  # Destination file exists
  if cmp -s ${srcPath} ${dstPath}; then
    # Destination file and source file are similar, exiting
    exit 0
  fi
fi

# Prepare destination directory
dstDir=$(dirname "${dstPath}")
mkdir -p ${dstDir}

# Copy file
cp -f ${srcPath} ${dstPath}.tmp

# Get file extension
srcFile=$(basename "${srcPath}")
srcExt="${srcFile##*.}"

if test "${srcExt}" = "h" -o "${srcExt}" = "cc"; then
  # Update content
  sed -i -e s/"oops::Variables"/"oops::patch::Variables"/g ${dstPath}.tmp
  sed -i -e s/" Variables"/" patch::Variables"/g ${dstPath}.tmp
  sed -i -e s/"<Variables>"/"<patch::Variables>"/g ${dstPath}.tmp
  sed -i -e s/"class patch::Variables;"/"namespace patch{\nclass Variables;\n}"/g ${dstPath}.tmp
  sed -i -e s/"^Variables"/"patch::Variables"/g ${dstPath}.tmp
fi

if test -f "${dstPath}.patch"; then
  # Apply residual patch
  patch -s ${dstPath}.tmp ${dstPath}.patch

  # Compare and update if needed
  if cmp -s ${dstPath}.tmp ${dstPath}; then
    rm -f ${dstPath}.tmp
    echo "--  - Update not needed for: ${dstPath}"
  else
    mv ${dstPath} ${dstPath}.bak
    mv ${dstPath}.tmp ${dstPath}
    echo "--  - Update needed for: ${dstPath}"
    echo "--      Diff command: diff -u ${dstPath} ${dstPath}.bak"
 fi
else
  if test -f "${dstPath}"; then
    # Create patch if needed
    diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch
    patchLength=`cat ${dstPath}.patch | wc -l`
    if test "${patchLength}" = "0"; then 
      rm -f ${dstPath}.patch
    else
      echo "--  - Patch created for: "${dstPath}
      echo "--      "`diffstat -C ${dstPath}.patch | head -n 1`
    fi
    rm -f ${dstPath}.tmp
  else
    # New file
    mv ${dstPath}.tmp ${dstPath}
  fi
fi
