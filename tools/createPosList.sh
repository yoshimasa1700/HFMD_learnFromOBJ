#!/bin/zsh

basename=$1
echo "${basename}"
classFolderList=(`echo ${basename}/*`)

for cFolder in ${classFolderList}
do
    if [ -d $cFolder ]; then
	#classes+=${cFolder##*/}
	echo ${cFolder##*/}

	instanceFolderList=(`echo ${cFolder}/*`)

	for iFolder in ${instanceFolderList}
	do
	    echo $iFolder
	    
	    fileList=(`echo ${iFolder}/*`)

	    rm ${iFolder}/dataList2.0.txt
	    rm ${iFolder}/dataListInst2.0.txt
	    rm ${iFolder}/testDataList2.0.txt
	    rm ${iFolder}/testDataListInst2.0.txt
	    touch ${iFolder}/dataList2.0.txt
	    touch ${iFolder}/dataListInst2.0.txt
	    touch ${iFolder}/testDataList2.0.txt
	    touch ${iFolder}/testDataListInst2.0.txt


	    num=0
	    for f in ${fileList}
	    do
		if [ ! -d $f ]; then
		    
		    case ${f##*/} in
			*_crop*)
			    bName=${f##*/}
			    bbName=${bName%_*}
			    #echo ${bbName}
			    angle=`expr ${bbName##*_} \* 6`
			    #echo ${angle}

			    width=`identify -format "%w" ${iFolder}/${bName%_*}_crop.png`
			    height=`identify -format "%h" ${iFolder}/${bName%_*}_crop.png`
			    #echo "$width"
			    centx=`expr $width / 2`
			    centy=`expr $height / 2`

			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png ${bName%_*}_maskcrop.png ${cFolder##*/} 0 0 $angle >> ${iFolder}/dataList2.0.txt
			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png ${bName%_*}_maskcrop.png ${iFolder##*/} 0 0 $angle >> ${iFolder}/dataListInst2.0.txt
			    echo "${bName%_*}_crop.png ${bName%_*}_depthcrop.png nodata ${cFolder##*/} ${centx} ${centy} 0 0 ${angle} EOL" >> ${iFolder}/testDataList2.0.txt
			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png nodata ${iFolder##*/} ${centx} ${centy} 0 0 ${angle} EOL >> ${iFolder}/testDataListInst2.0.txt			
	#echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png nodata ${cFolder##*/} ${centx} ${centy} ${angle} EOL
			    #echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png nodata ${iFolder##*/} ${centx} ${centy} ${angle} EOL 		
		
num=`expr ${num} + 1`	
			
			
			
		esac

		#	echo ${num}
		fi
	    done
	#echo ${num}

	    #cat ${iFolder}/dataList2.0.txt | sed -i "1 ${num}"  > ${iFolder}/dataList2.0.txt
	    #sed -i "1s \\${num}" '${iFolder}/dataList2.0.txt'

	    sed -i"" "1s/^/${num}\n/" ${iFolder}/dataList2.0.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/dataListInst2.0.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/testDataList2.0.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/testDataListInst2.0.txt

	done
    fi
done

# if [ $# -ge 1 ]; then
#     echo "60" > dataList.txt
#     for i in `seq 1 60`
#     do
# 	numberedname=${babasename}_1_${i}_
# 	echo ${numberedname}crop.png ${numberedname}depthcrop.png >> negDataList.txt
# 	mv ./negDataList.txt ./dataset/${basename}/${babasename}
#     done
# else
#     echo "usage : ./createNegList.sh [basename *example bottle_1]"
#     exit 1
# fi
