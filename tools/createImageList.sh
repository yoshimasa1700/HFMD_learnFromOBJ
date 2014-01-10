#!/bin/zsh

basename=$1
echo "${basename}"
classFolderList=(`echo ${basename}/*`)

for cFolder in ${classFolderList}
do
    if [ -d $cFolder ]; then
	echo ${cFolder##*/}

	instanceFolderList=(`echo ${cFolder}/*`)

	for iFolder in ${instanceFolderList}
	do
	    echo $iFolder
	    
	    fileList=(`echo ${iFolder}/*`)

	    ## delete exsist txt file
	    rm -f ${iFolder}/trainImageList.txt
	    rm -f ${iFolder}/trainImageListI.txt

	    rm -f ${iFolder}/testImageList.txt
	    rm -f ${iFolder}/testImageList.txt

	    rm -f ${iFolder}/negImageList.txt

	    ## create new image list
	    touch ${iFolder}/trainImageList.txt
	    touch ${iFolder}/trainImageListI.txt

	    touch ${iFolder}/testImageList.txt
	    touch ${iFolder}/testImageList.txt
	    
	    touch ${iFolder}/negImageList.txt

	    num=0
	    for f in ${fileList}
	    do
		if [ ! -d $f ]; then
		    
		    case ${f##*/} in
			*_crop*)
			    bName=${f##*/}
			    bbName=${bName%_*}

			    angle=`expr ${bbName##*_} \* 6`

			    width=`identify -format "%w" ${iFolder}/${bName%_*}_crop.png`
			    height=`identify -format "%h" ${iFolder}/${bName%_*}_crop.png`

			    centx=`expr $width / 2`
			    centy=`expr $height / 2`

			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png NULL ${cFolder##*/} 0 0 $angle >> ${iFolder}/trainImageList.txt
			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png NULL ${iFolder##*/} 0 0 $angle >> ${iFolder}/trainImageListI.txt
			    echo "${bName%_*}_crop.png ${bName%_*}_depthcrop.png NULL ${cFolder##*/} ${centx} ${centy} 0 0 ${angle} EOL" >> ${iFolder}/testImageList.txt
			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png NULL ${iFolder##*/} ${centx} ${centy} 0 0 ${angle} EOL >> ${iFolder}/testImageListI.txt
			    echo ${bName%_*}_crop.png ${bName%_*}_depthcrop.png >> ${iFolder}/negImageList.txt
		
			    num=`expr ${num} + 1`
		esac
		fi
	    done

	    sed -i"" "1s/^/${num}\n/" ${iFolder}/trainImageList.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/trainImageListI.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/testImageList.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/testImageListI.txt
	    sed -i"" "1s/^/${num}\n/" ${iFolder}/negImageList.txt

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
