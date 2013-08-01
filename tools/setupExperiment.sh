#!/bin/zsh
echo "test"
echo $#
currentD=`pwd`
echo $currentD

needFile+=("learning")
needFile+=("objectPoseEstimation")
needFile+=("config.xml")
#needFile+=("changeConfig.sh")

needFolder+=("dataset")
needFolder+=("negdata")
needFolder+=("trees")
needFolder+=("testdata")
needFolder+=("tools")
needFolder+=("models")

echo ${needFile[@]}
echo ${needFolder[@]}

if [ $# -ge 1 ]; then
	echo "File name is $1"
	cd ../
	mkdir $1
	cd $1
	for filename in "${needFile[@]}"
	do
		filenamepath=${currentD}/${filename}
		cp $filenamepath ./
	done
	for foldername in "${needFolder[@]}"
	do
		folderpath=${currentD}/${foldername}
		echo ${folderpath}
		cp -r $folderpath ./
	done
else
	echo "argument 1 needed"
	exit 1
fi

