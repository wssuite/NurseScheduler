#!/bin/bash

# generate script
source generateScript.sh $1 $2

# run script
echo "Run: ./${scriptfile}"
chmod +x *.jar
cp validator.jar ./bin
./${scriptfile}

# if $3 defined, test the total cost
if [ -z $3 ]
then
	exit 0
fi

# fetch the total cost and check if is the right one
cat ${outputDir}/validator.txt
result=$(cat ${outputDir}/validator.txt | grep "Total cost: $3")
if [ -z "$result" ]
then
  echo "error"
  exit 1
fi
exit 0
