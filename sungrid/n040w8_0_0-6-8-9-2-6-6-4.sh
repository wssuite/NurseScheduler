#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n040w8 0 0 6 8 9 2 6 6 4 n040w8_0_0-6-8-9-2-6-6-4 > outfiles/Competition/n040w8_0_0-6-8-9-2-6-6-4/log.txt

instance=n040w8
weeksValue=(0 6 8 9 2 6 6 4 )

demand0="WD-${instance}-"
solutionFile="outfiles/n040w8_0_0-6-8-9-2-6-6-4/OptSol-n040w8-06892664-"
weeks=""
sols=""
i=0

for var in ${weeksValue[*]}
do
demand[$i]="datatsets/${instance}/${demand0}${var}.txt"
weeks="${weeks} ${demand[$i]}"
solution[$i]="${solutionFile}${var}-${i}.txt"
sols="${sols} ${solution[$i]}"
((i++))
done

java -jar validator.jar --sce datasets/n040w8/Sc-n040w8.txt --his datasets/n040w8/H0-n040w8-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n040w8_0_0-6-8-9-2-6-6-4/validatorOutput.txt 

exit 0;
