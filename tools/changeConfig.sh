#!/bin/zsh

alias xml='xmlstarlet'
if [ $# -ge 3 ];then
	echo "config file name is $1"
	echo "attribute is $2"
	echo "value is $3"
	echo "xpath is root/"

	xml ed -u "/root/$2" -v $3 $1 > tempConfig.xml
	cp tempConfig.xml $1
	rm tempConfig.xml
else
	echo "usage: changeConfig [configfile] [xpath] [value]"
	exit 1
fi
