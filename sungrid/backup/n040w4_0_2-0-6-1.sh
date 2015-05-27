#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n040w4 0 2 0 6 1 n040w4_0_2-0-6-1 $1 > outfiles/Competition/n040w4_0_2-0-6-1/${1}Log.txt

instance=n040w4
weeksValue=(2 0 6 1 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n040w4_0_2-0-6-1/${1}Sol-n040w4-2061-"
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

java -jar validator.jar --sce datasets/n040w4/Sc-n040w4.txt --his datasets/n040w4/H0-n040w4-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n040w4_0_2-0-6-1/${1}ValidatorOutput.txt 

exit 0;
