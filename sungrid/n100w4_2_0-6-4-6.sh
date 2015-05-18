#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n100w4 2 0 6 4 6 n100w4_2_0-6-4-6 > outfiles/Competition/n100w4_2_0-6-4-6/log.txt

instance=n100w4
weeksValue=(0 6 4 6 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n100w4_2_0-6-4-6/OptSol-n100w4-0646-"
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

java -jar validator.jar --sce datasets/n100w4/Sc-n100w4.txt --his datasets/n100w4/H0-n100w4-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n100w4_2_0-6-4-6/validatorOutput.txt 

exit 0;
