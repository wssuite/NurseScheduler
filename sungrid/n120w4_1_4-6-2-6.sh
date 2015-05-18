#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n120w4 1 4 6 2 6 n120w4_1_4-6-2-6 > outfiles/Competition/n120w4_1_4-6-2-6/log.txt

instance=n120w4
weeksValue=(4 6 2 6 )

demand0="WD-${instance}-"
solutionFile="outfiles/n120w4_1_4-6-2-6/OptSol-n120w4-4626-"
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

java -jar validator.jar --sce datasets/n120w4/Sc-n120w4.txt --his datasets/n120w4/H0-n120w4-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n120w4_1_4-6-2-6/validatorOutput.txt 

exit 0;
