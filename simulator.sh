#!/bin/bash
# simulator script: launch the simulator

dataFile="testdatasets/"
instance="n005w4"
solutionFile="Solution_H_0-WD_1-2-3-3/"
weeksValue=(1 2 3 3)

history="${dataFile}${instance}/H0-${instance}-0.txt"
scenario="${dataFile}${instance}/Sc-${instance}.txt"

demand0="WD-${instance}-"
weeks=""
i=0

for var in ${weeksValue[*]}
do
demand[$i]="${dataFile}${instance}/${demand0}${var}.txt"
((i++))
done

java -jar Simulator.jar --his $history --sce $scenario --weeks $weeks --solver ./roster --runDir bin/ --outDir outfiles/ 

exit 0;
