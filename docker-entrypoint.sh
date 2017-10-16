#!/bin/bash

# parse inputs
# parse the input competition instance name
echo "Instance: $1"
parse=(`echo $1 | tr '_' ' ' `)
instance=${parse[0]}
hist=${parse[1]}
weeks=${parse[2]}

# create the root of output directory if it does not exist
printf -v currenttime '%(%s)T' -1
outputDir="outfiles/$1/$currenttime"
echo "Create output directory: ${outputDir}"
mkdir -p "${outputDir}"

# set default timeout
timeout=$TIMEOUT
if [ -z $3 ]
then
	if [ -z "$TIMEOUT" ]
	then
	  timeout=30
	fi
else
	timeout=$3
fi

# set default param file
param=$PARAM
if [ -z $4 ]
then
	if [ -z "$PARAM" ]
	then
	  param="default.txt"
	fi
else
	param=$4
fi

# run the scheduler
echo "Run: ./bin/staticscheduler --dir datasets/ --instance $instance --weeks $weeks --his $hist --param paramfiles/$param --sol $outputDir --timeout $timeout"
./bin/staticscheduler --dir datasets/ --instance $instance --weeks $weeks --his $hist --param paramfiles/$param --sol $outputDir --timeout $timeout

# run the validator
./validator.sh $instance $weeks $hist $outputDir

# if $2 defined, test the total cost
if [ -z $2 ]
then
	exit 0
fi

# fetch the total cost and check if is the right one
cat ${outputDir}/validator.txt
result=$(cat ${outputDir}/validator.txt | grep "Total cost: $2")
if [ -z "$result" ]
then
  echo "error"
  exit 1
fi
exit 0
