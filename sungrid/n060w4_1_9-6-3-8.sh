#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w4 1 9 6 3 8 n060w4_1_9-6-3-8 > outfiles/Competition/n060w4_1_9-6-3-8/log.txt

instance=n060w4
weeksValue=(9 6 3 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/n060w4_1_9-6-3-8/OptSol-n060w4-9638-"
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

java -jar validator.jar --sce datasets/n060w4/Sc-n060w4.txt --his datasets/n060w4/H0-n060w4-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n060w4_1_9-6-3-8/validatorOutput.txt 

exit 0;
