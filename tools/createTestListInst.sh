#!/bin/zsh

basename=$1
number=$2
echo "${basename}_${number}"
babasename=${basename}_${number}
bababasename=${basename}${number}

if [ $# -ge 1 ]; then
    echo "60" > testDataListInst.txt
    for i in `seq 1 60`
    do
	numberedname=${babasename}_1_${i}_
	width=`identify -format "%w" ./dataset/${basename}/${babasename}/${numberedname}crop.png`
	height=`identify -format "%h" ./dataset/${basename}/${babasename}/${numberedname}crop.png`
	#echo "$width"
	centx=`expr $width / 2`
	centy=`expr $height / 2`
	angle=`expr $i \* 6`

	echo ${numberedname}crop.png ${numberedname}depthcrop.png nodata ${bababasename} ${centx} ${centy} ${angle} EOL >> testDataListInst.txt
	
    done
    mv ./testDataListInst.txt ./dataset/${basename}/${babasename}
else
    echo "usage : ./createNegList.sh [basename *example bottle_1]"
    exit 1
fi
