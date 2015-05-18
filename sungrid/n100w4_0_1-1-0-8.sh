#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n100w4 0 1 1 0 8 n100w4_0_1-1-0-8 > outfiles/Competition/n100w4_0_1-1-0-8/log.txt

instance=n100w4
weeksValue=(1 1 0 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/n100w4_0_1-1-0-8/OptSol-n100w4-1108-"
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

java -jar validator.jar --sce datasets/n100w4/Sc-n100w4.txt --his datasets/n100w4/H0-n100w4-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n100w4_0_1-1-0-8/validatorOutput.txt 

exit 0;
