#!/bin/sh
filename=$1
cat ${filename} | while read line
do
  ./tools/createTestListInst.sh ${line}
done 
