#!/bin/bash
# validator script: launch the validator

dataFile="datasets/"
instance=$1
solutionDir="outfiles/"
weeksName=$2
historyIndex=$3
history="${dataFile}${instance}/H0-${instance}-${historyIndex}.txt"
scenario="${dataFile}${instance}/Sc-${instance}.txt"
solutionFile="outfiles/$4/"

demand0="WD-${instance}-"
solution0="sol-week"
weeks=""
sols=""
i=0

for var in ${weeksName//-/ }
do
demand[$i]="${dataFile}${instance}/${demand0}${var}.txt"
weeks="${weeks} ${demand[$i]}"

solution[$i]="${solutionFile}${solution0}${i}.txt"
sols="${sols} ${solution[$i]}"

((i++))
done

echo "java -jar validator.jar --his $history --sce $scenario --weeks $weeks --sols $sols"
val=$(java -jar validator.jar --his $history --sce $scenario --weeks $weeks --sols $sols)
echo $val
echo $val >> "${solutionFile}validator.txt"
exit 0;