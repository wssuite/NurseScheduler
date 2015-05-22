#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w4 1 6 1 1 5 n060w4_1_6-1-1-5 $1 > outfiles/Competition/n060w4_1_6-1-1-5/${1}Log.txt

instance=n060w4
weeksValue=(6 1 1 5 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n060w4_1_6-1-1-5/${1}Sol-n060w4-6115-"
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

java -jar validator.jar --sce datasets/n060w4/Sc-n060w4.txt --his datasets/n060w4/H0-n060w4-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n060w4_1_6-1-1-5/${1}ValidatorOutput.txt 

exit 0;
