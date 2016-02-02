#!/bin/bash

cd `dirname $0`

export MYTRIM_SEED=39172
../../build/mytrim_uo2 out 10 0.1 1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
  echo "FAILED to run mytrim_uo2!"
  exit 1
fi

NFAIL=0
for OUT in out.Erec out.clcoor out.dist
do
  ../csv_diff.sh $OUT 0.0001
  if [ $? -ne 0 ]
  then
    echo "Difference in $OUT"
    let NFAIL=NFAIL+1
  fi
done

if [ $NFAIL -gt 0 ]; then
  echo "FAILED!"
  exit 1
fi
