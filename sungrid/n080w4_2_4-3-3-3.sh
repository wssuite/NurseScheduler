#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n080w4 2 4 3 3 3 n080w4_2_4-3-3-3 $1 > outfiles/Competition/n080w4_2_4-3-3-3/${1}Log.txt

instance=n080w4
weeksValue=(4 3 3 3 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n080w4_2_4-3-3-3/${1}Sol-n080w4-4333-"
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

java -jar validator.jar --sce datasets/n080w4/Sc-n080w4.txt --his datasets/n080w4/H0-n080w4-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n080w4_2_4-3-3-3/${1}ValidatorOutput.txt 

exit 0;
