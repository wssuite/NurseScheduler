#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n100w8 0 0 1 7 8 9 1 5 4 n100w8_0_0-1-7-8-9-1-5-4 > outfiles/Competition/n100w8_0_0-1-7-8-9-1-5-4/log.txt

instance=n100w8
weeksValue=(0 1 7 8 9 1 5 4 )

demand0="WD-${instance}-"
solutionFile="outfiles/n100w8_0_0-1-7-8-9-1-5-4/OptSol-n100w8-01789154-"
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

java -jar validator.jar --sce datasets/n100w8/Sc-n100w8.txt --his datasets/n100w8/H0-n100w8-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n100w8_0_0-1-7-8-9-1-5-4/validatorOutput.txt 

exit 0;
