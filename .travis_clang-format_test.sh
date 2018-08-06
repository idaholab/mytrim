#!/bin/bash

# apply clang-format
find . -name \*.C -o -name \*.h -exec clang-format -i \{\} \;

# see if any changes were made
git diff | tee .diff | grep diff > /dev/null || exit 0

# changes were made -> branch is not formatted according to rules
cat .diff
exit 1
