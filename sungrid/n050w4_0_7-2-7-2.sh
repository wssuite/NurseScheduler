#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n050w4 0 7 2 7 2 n050w4_0_7-2-7-2 > outfiles/Competition/n050w4_0_7-2-7-2/log.txt

instance=n050w4
weeksValue=(7 2 7 2 )

demand0="WD-${instance}-"
solutionFile="outfiles/n050w4_0_7-2-7-2/OptSol-n050w4-7272-"
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

java -jar validator.jar --sce datasets/n050w4/Sc-n050w4.txt --his datasets/n050w4/H0-n050w4-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n050w4_0_7-2-7-2/validatorOutput.txt 

exit 0;
