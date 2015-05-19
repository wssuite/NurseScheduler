#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n100w8 1 2 4 7 9 3 9 2 8 n100w8_1_2-4-7-9-3-9-2-8 > outfiles/Competition/n100w8_1_2-4-7-9-3-9-2-8/log.txt

instance=n100w8
weeksValue=(2 4 7 9 3 9 2 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n100w8_1_2-4-7-9-3-9-2-8/OptSol-n100w8-24793928-"
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

java -jar validator.jar --sce datasets/n100w8/Sc-n100w8.txt --his datasets/n100w8/H0-n100w8-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n100w8_1_2-4-7-9-3-9-2-8/validatorOutput.txt 

exit 0;
