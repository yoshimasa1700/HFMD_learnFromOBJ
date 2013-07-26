#!/bin/zsh

if [ $# -ge 5 ]; then
basename=$1
param=$2
min=$3
max=$4
stride=$5

echo "basename is $1"
echo "param is $2"
echo "min is $3 max is $4 stride is $5"

sa=`expr $4 - $3`
haba=`expr $sa / $5`
currentDir= pwd
echo "$currentDir"
echo "sa is $haba"
	for i in `seq 0 $haba`
	do
		echo "$i"
		echo "$currentDir"
		lerolero=`expr $i \* $stride`
		lerolero2=`expr $lerolero + $min`
		foldername=$basename$param$lerolero2
		echo "$lerolero2"
		echo "$foldername"
		#./tools/setupExperiment.sh $foldername
#cp ./tools/analyzeResult.py ../${foldername}/tools
#cd ../$foldername
		#./tools/changeConfig.sh config.xml $param $lerolero2
		#./learning
		#./tools/changeConfig.sh config.xml stride 5
		#./objectPoseEstimation
        ./analyzeResult.py ./${foldername}Res.txt ./${foldername}Ana.txt

		cd "$currentDir" 
		echo "$currentDir" 
	done
else
	echo "usage : paramChangeExam.sh [basename] [param] [min] [max] [stride]"
	exit 1
fi
