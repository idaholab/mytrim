#!/bin/bash

FILE=$1
TOL=$2

# detect number of colums
NCOL=`awk '{ print NF; exit}' $FILE`

# paste and compare each column
for COL in `seq $NCOL`; do
  let COL2=COL+NCOL
  OUT=`paste $FILE gold/$FILE | awk '{ d=$'$COL'-$'$COL2'; if (d > '$TOL' || d < -'$TOL') { print "FAIL"; exit; } }'`
  if [ x$OUT == "xFAIL" ]; then
    exit 1
  fi
done
