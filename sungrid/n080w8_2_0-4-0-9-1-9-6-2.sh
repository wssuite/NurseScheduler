#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n080w8 2 0 4 0 9 1 9 6 2 n080w8_2_0-4-0-9-1-9-6-2 > outfiles/Competition/n080w8_2_0-4-0-9-1-9-6-2/log.txt

instance=n080w8
weeksValue=(0 4 0 9 1 9 6 2 )

demand0="WD-${instance}-"
solutionFile="outfiles/n080w8_2_0-4-0-9-1-9-6-2/OptSol-n080w8-04091962-"
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

java -jar validator.jar --sce datasets/n080w8/Sc-n080w8.txt --his datasets/n080w8/H0-n080w8-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n080w8_2_0-4-0-9-1-9-6-2/validatorOutput.txt 

exit 0;
