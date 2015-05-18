#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n080w4 2 6 0 4 8 n080w4_2_6-0-4-8 > outfiles/Competition/n080w4_2_6-0-4-8/log.txt

instance=n080w4
weeksValue=(6 0 4 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/n080w4_2_6-0-4-8/OptSol-n080w4-6048-"
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

java -jar validator.jar --sce datasets/n080w4/Sc-n080w4.txt --his datasets/n080w4/H0-n080w4-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n080w4_2_6-0-4-8/validatorOutput.txt 

exit 0;
