#!/bin/zsh

basename=$1
number=$2
echo "${basename}_${number}"
babasename=${basename}_${number}

if [ $# -ge 1 ]; then
    echo "60" > negDataList.txt
    for i in `seq 1 60`
    do
	numberedname=${babasename}_1_${i}_
	echo ${numberedname}crop.png ${numberedname}depthcrop.png >> negDataList.txt
	mv ./negDataList.txt ./dataset/${basename}/${babasename}
    done
else
    echo "usage : ./createNegList.sh [basename *example bottle_1]"
    exit 1
fi
