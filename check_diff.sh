#!/usr/bin/env bash

# Source and destination paths
srcPath=$1
dstPath=$2

# Make destination directory
dstDir=$(dirname "${dstPath}")
mkdir -p ${dstDir}

# Get file extension
srcFile=$(basename "${srcPath}")
srcExt="${srcFile##*.}"

# Copy file
cp -f ${srcPath} ${dstPath}.tmp

# Update content
sed -i -e s/"oops::Variables"/"oops::patch::Variables"/g ${dstPath}.tmp
sed -i -e s/" Variables"/" patch::Variables"/g ${dstPath}.tmp
sed -i -e s/"<Variables>"/"<patch::Variables>"/g ${dstPath}.tmp
sed -i -e s/"class patch::Variables;"/"namespace patch{\nclass Variables;\n}"/g ${dstPath}.tmp
sed -i -e s/"^Variables"/"patch::Variables"/g ${dstPath}.tmp

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
  # Create patch for new file
  diff -u ${dstPath}.tmp ${dstPath} > ${dstPath}.patch
  echo "--  - Patch created for: ${dstPath}"
fi
