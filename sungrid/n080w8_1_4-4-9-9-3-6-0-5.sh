#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n080w8 1 4 4 9 9 3 6 0 5 n080w8_1_4-4-9-9-3-6-0-5 > outfiles/Competition/n080w8_1_4-4-9-9-3-6-0-5/log.txt

instance=n080w8
weeksValue=(4 4 9 9 3 6 0 5 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n080w8_1_4-4-9-9-3-6-0-5/OptSol-n080w8-44993605-"
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

java -jar validator.jar --sce datasets/n080w8/Sc-n080w8.txt --his datasets/n080w8/H0-n080w8-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n080w8_1_4-4-9-9-3-6-0-5/validatorOutput.txt 

exit 0;
