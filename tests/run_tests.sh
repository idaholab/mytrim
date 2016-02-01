#!/bin/bash

export MYTRIM_DATADIR=$TRAVIS_BUILD_DIR/data

for TEST in `find . -name test.sh`
do
  echo running $TEST ...
  $TEST
  if [ $? -ne 0 ]
  then
    exit 1
  fi
done

echo "ALL PASSED."
