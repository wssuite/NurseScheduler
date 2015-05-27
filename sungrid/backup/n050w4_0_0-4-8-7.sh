#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n050w4 0 0 4 8 7 n050w4_0_0-4-8-7 $1 > outfiles/Competition/n050w4_0_0-4-8-7/${1}Log.txt

instance=n050w4
weeksValue=(0 4 8 7 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n050w4_0_0-4-8-7/${1}Sol-n050w4-0487-"
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

java -jar validator.jar --sce datasets/n050w4/Sc-n050w4.txt --his datasets/n050w4/H0-n050w4-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n050w4_0_0-4-8-7/${1}ValidatorOutput.txt 

exit 0;
