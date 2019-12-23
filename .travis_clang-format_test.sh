#!/bin/bash

CLANG_FORMAT=clang-format
clang-format-9 --version >/dev/null 2>&1 && CLANG_FORMAT=clang-format-9

# output version
$CLANG_FORMAT --version

# apply clang-format
find . -name \*.C -o -name \*.h -exec $CLANG_FORMAT -i \{\} \;

# see if any changes were made
git diff | tee .diff | grep diff > /dev/null
if [ $? != 0 ]
then
  echo Passed.
  exit 0
fi

# changes were made -> branch is not formatted according to rules
cat .diff
exit 1
