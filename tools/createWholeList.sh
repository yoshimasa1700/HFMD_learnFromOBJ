#!/bin/sh
filename=$1
cat ${filename} | while read line
do
  ./tools/createTestList.sh ${line}
done 
