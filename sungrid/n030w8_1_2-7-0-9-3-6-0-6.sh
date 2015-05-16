#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n030w8 1 2 7 0 9 3 6 0 6 n030w8_1_2-7-0-9-3-6-0-6 > outfiles/Competition/n030w8_1_2-7-0-9-3-6-0-6/log.txt

instance=n030w8
weeksValue=(2 7 0 9 3 6 0 6 )

demand0="WD-${instance}-"
solutionFile="outfiles/n030w8_1_2-7-0-9-3-6-0-6/OptSol-n030w8-27093606-"
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

java -jar validator.jar --sce datasets/n030w8/Sc-n030w8.txt --his datasets/n030w8/H0-n030w8-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n030w8_1_2-7-0-9-3-6-0-6/validatorOutput.txt 

exit 0;
