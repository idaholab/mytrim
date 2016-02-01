#!/bin/bash

cd `dirname $0`
ls -la ../../build

export MYTRIM_SEED=39172
../../build/mytrim_uo2 out 10 0.1 1
#> /dev/null 2>&1

for OUT in out.Erec out.clcoor out.dist
do
  diff $OUT gold/$OUT
  #> /dev/null
  if [ $? -ne 0 ]
  then
    echo "FAILED! difference in $OUT"
    exit 1
  fi
done
