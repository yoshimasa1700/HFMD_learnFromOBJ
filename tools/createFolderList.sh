#!/bin/zsh

basename=$1
echo "${basename}"
classFolderList=(`echo ${basename}/*`)

rm ${basename}/trainDataList2.0.txt
rm ${basename}/testDataList2.0.txt
touch ${basename}/trainDataList2.0.txt
touch ${basename}/testDataList2.0.txt

num=0

for cFolder in ${classFolderList}
do
    if [ -d $cFolder ]; then
	#classes+=${cFolder##*/}
	echo ${cFolder##*/}

	instanceFolderList=(`echo ${cFolder}/*`)

	for iFolder in ${instanceFolderList}
	do
	    echo $iFolder
	    
	    short=${iFolder#*/}
	    echo ${short#*/} >> ${basename}/trainDataList2.0.txt
	    echo ${short#*/} >> ${basename}/testDataList2.0.txt
	    #cat ${iFolder}/dataList2.0.txt | sed -i "1 ${num}"  > ${iFolder}/dataList2.0.txt
	    #sed -i "1s \\${num}" '${iFolder}/dataList2.0.txt'
	    num=`expr $num + 1`
	done
    fi
    echo $num
        
done

sed -i "1s/^/${num}\n/" ${basename}/trainDataList2.0.txt
    sed -i "1s/^/${num}\n/" ${basename}/testDataList2.0.txt

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
