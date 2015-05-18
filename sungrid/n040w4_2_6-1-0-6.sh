#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n040w4 2 6 1 0 6 n040w4_2_6-1-0-6 > outfiles/Competition/n040w4_2_6-1-0-6/log.txt

instance=n040w4
weeksValue=(6 1 0 6 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n040w4_2_6-1-0-6/OptSol-n040w4-6106-"
weeks=""
sols=""
i=0

for var in ${weeksValue[*]}
do
demand[$i]="datasets/${instance}/${demand0}${var}.txt"
weeks="${weeks} ${demand[$i]}"
solution[$i]="${solutionFile}${var}-${i}.txt"
sols="${sols} ${solution[$i]}"
((i++))
done

java -jar validator.jar --sce datasets/n040w4/Sc-n040w4.txt --his datasets/n040w4/H0-n040w4-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n040w4_2_6-1-0-6/validatorOutput.txt 

exit 0;
