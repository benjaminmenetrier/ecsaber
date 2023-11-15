#!/usr/bin/env bash

# Source and destination paths
srcPath=$1
dstPath=$2
commandPath=$3

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
  cp -f ${dstPath}.tmp ${dstPath}.tmp.bak
  patch -s ${dstPath}.tmp ${dstPath}.patch

  # Compare and update if needed
  if cmp -s ${dstPath}.tmp ${dstPath}; then
    rm -f ${dstPath}.tmp ${dstPath}.tmp.bak
  else
    echo "--  - Update needed for: ${dstPath}"
    echo "meld ${dstPath}.tmp ${dstPath} ${dstPath}.tmp.bak; diff -u ${dstPath}.tmp.bak ${dstPath} > ${dstPath}.patch; rm -f ${dstPath}.tmp ${dstPath}.tmp.bak" >> ${commandPath}
    echo "" >> ${commandPath}
 fi
else
  if test -f "${dstPath}"; then
    # Create patch if needed
    diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch
    patchLength=`cat ${dstPath}.patch | wc -l`
    if test "${patchLength}" = "0"; then 
      rm -f ${dstPath}.patch
    else
      echo "--  - New patch needed for: "${dstPath}
      rm -f ${dstPath}.patch
      echo "meld ${dstPath}.tmp ${dstPath}; diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch; patchLength=\`cat ${dstPath}.patch | wc -l\`; if test "${patchLength}" = "0"; then rm -f ${dstPath}.patch;fi;rm -f ${dstPath}.tmp" >> ${commandPath}
      echo "" >> ${commandPath}
    fi
  else
    # New file
    mv ${dstPath}.tmp ${dstPath}
  fi
fi
