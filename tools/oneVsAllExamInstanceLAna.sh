#!/bin/zsh

classfilename=$1
folderfilename=$2

cat ${classfilename}| while read line
do
	classname+=(${line})
done

cat ${folderfilename}| while read line
do
    foldername+=(`echo "${line}"`)
done

for i in `seq 1 ${#classname[@]}`
do
    #実験フォルダの準備
    #./tools/changeConfig.sh ./config.xml stride 1
    #./tools/setupExperiment.sh ${classname[$i]}VsAll
    cd ../${classname[$i]}VsAll
    ./tools/changeConfig.sh ./config.xml stride 1
    echo ${classname[$i]}
    # num=`expr ${#foldername[@]} / 2`
    
    # posfolder=0
    # negfolder=0
    
    # pwd
    
    # for j in `seq 1 $num`
    # do
	
    # 	k=`expr $j\*2`
    # 	if [ "${classname[$i]}" = "${foldername[$k-1]}_${foldername[$k]}" ]; then
    # 	    echo ${foldername[$k-1]}_${foldername[$k]}
    # 	    posfolder=`expr $posfolder + 1`
    # 	fi
    # 	if [ "${classname[$i]}" != "${foldername[$k-1]}_${foldername[$k]}" ]; then
    # 	    echo ${foldername[$k-1]}_${foldername[$k]}
    # 	    negfolder=`expr $negfolder + 1`
    # 	fi
    # done
    
    # echo ${posfolder} > ./trainData.txt
    # echo ${negfolder} > ./negDataFolderList.txt
    # for j in `seq 1 $num`
    # do
	
    # 	k=`expr $j\*2`
    # 	if [ "${classname[$i]}" = "${foldername[$k-1]}_${foldername[$k]}" ]; then
    # 	    echo "${foldername[$k-1]}/${foldername[$k-1]}_${foldername[$k]}" >> trainData.txt
	    
    # 	fi
    # 	if [ "${classname[$i]}" != "${foldername[$k-1]}_${foldername[$k]}" ]; then
    # 	    echo "../dataset/${foldername[$k-1]}/${foldername[$k-1]}_${foldername[$k]}" >> negDataFolderList.txt
	    
    # 	fi
    # done
    # 	#pwd
    # 	#ls
    
    #mv ./trainData.txt ./dataset
    #mv ./negDataFolderList.txt ./negdata
    
    expDir="../${classname[$i]}VsAll"

pwd
#./learning 
./tools/changeConfig.sh ./config.xml stride 5
#./objectPoseEstimation
cp -r ../base/dataset ./
cp ../base/testdata/testData.txt ./testdata
cp ../base/config.xml ./
../base/tools/analyzeResult.py ./detectionResult.txt ../base/${classname[$i]}Ana.txt ../base/${classname[$i]}Res.txt
cd ../base
done


