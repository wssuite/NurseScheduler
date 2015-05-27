#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w8 2 1 0 3 4 0 3 9 1 n060w8_2_1-0-3-4-0-3-9-1 $1 > outfiles/Competition/n060w8_2_1-0-3-4-0-3-9-1/${1}Log.txt

instance=n060w8
weeksValue=(1 0 3 4 0 3 9 1 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n060w8_2_1-0-3-4-0-3-9-1/${1}Sol-n060w8-10340391-"
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

java -jar validator.jar --sce datasets/n060w8/Sc-n060w8.txt --his datasets/n060w8/H0-n060w8-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n060w8_2_1-0-3-4-0-3-9-1/${1}ValidatorOutput.txt 

exit 0;
