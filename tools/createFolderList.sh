#!/bin/zsh

basename=$1
echo "${basename}"
classFolderList=(`echo ${basename}/*`)

rm ${basename}/trainFolderList.txt
rm ${basename}/testFolderList.txt
touch ${basename}/trainFolderList.txt
touch ${basename}/testFolderList.txt

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
	    echo ${short#*/} >> ${basename}/trainFolderList.txt
	    echo ${short#*/} >> ${basename}/testFolderList.txt
	    #cat ${iFolder}/dataList2.0.txt | sed -i "1 ${num}"  > ${iFolder}/dataList2.0.txt
	    #sed -i "1s \\${num}" '${iFolder}/dataList2.0.txt'
	    num=`expr $num + 1`
	done
    fi
    echo $num
        
done

sed -i "1s/^/${num}\n/" ${basename}/trainFolderList.txt
    sed -i "1s/^/${num}\n/" ${basename}/testFolderList.txt

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
