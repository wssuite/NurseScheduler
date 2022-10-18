#!/bin/bash
# validator script: launch the validator

dataFile="datasets/INRC2/"
instance=$1
solutionDir="outfiles/"
weeksName=$2
historyIndex=$3
history="${dataFile}${instance}/H0-${instance}-${historyIndex}.txt"
scenario="${dataFile}${instance}/Sc-${instance}.txt"

demand0="WD-${instance}-"
solution0="sol-week"
weeks=""
sols=""
i=0

for var in ${weeksName//-/ }
do
demand[$i]="${dataFile}${instance}/${demand0}${var}.txt"
weeks="${weeks} ${demand[$i]}"

solution[$i]="$4/${solution0}${i}.txt"
sols="${sols} ${solution[$i]}"

((i++))
done

echo "java -jar validator.jar --his $history --sce $scenario --weeks $weeks --sols $sols"
java -jar validator.jar --his $history --sce $scenario --weeks $weeks --sols $sols > "$4/validator.txt"
exit 0;
