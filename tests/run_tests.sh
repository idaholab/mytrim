#!/bin/bash

if [ ! -z $TRAVIS_BUILD_DIR ]; then
  export MYTRIM_DATADIR=$TRAVIS_BUILD_DIR/data
fi

for TEST in `find . -name test.sh`; do
  echo running $TEST ...
  $TEST
  if [ $? -ne 0 ]; then
    exit 1
  fi
done

echo "ALL PASSED."
