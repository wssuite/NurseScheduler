#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n050w8 1 9 7 5 3 8 8 3 1 n050w8_1_9-7-5-3-8-8-3-1 $1 > outfiles/Competition/n050w8_1_9-7-5-3-8-8-3-1/${1}Log.txt

instance=n050w8
weeksValue=(9 7 5 3 8 8 3 1 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n050w8_1_9-7-5-3-8-8-3-1/${1}Sol-n050w8-97538831-"
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

java -jar validator.jar --sce datasets/n050w8/Sc-n050w8.txt --his datasets/n050w8/H0-n050w8-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n050w8_1_9-7-5-3-8-8-3-1/${1}ValidatorOutput.txt 

exit 0;
