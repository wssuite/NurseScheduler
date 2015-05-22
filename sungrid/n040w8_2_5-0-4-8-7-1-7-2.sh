#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n040w8 2 5 0 4 8 7 1 7 2 n040w8_2_5-0-4-8-7-1-7-2 $1 > outfiles/Competition/n040w8_2_5-0-4-8-7-1-7-2/${1}Log.txt

instance=n040w8
weeksValue=(5 0 4 8 7 1 7 2 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n040w8_2_5-0-4-8-7-1-7-2/${1}Sol-n040w8-50487172-"
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

java -jar validator.jar --sce datasets/n040w8/Sc-n040w8.txt --his datasets/n040w8/H0-n040w8-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n040w8_2_5-0-4-8-7-1-7-2/${1}ValidatorOutput.txt 

exit 0;
