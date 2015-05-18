#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w8 0 6 2 9 9 0 8 1 3 n060w8_0_6-2-9-9-0-8-1-3 > outfiles/Competition/n060w8_0_6-2-9-9-0-8-1-3/log.txt

instance=n060w8
weeksValue=(6 2 9 9 0 8 1 3 )

demand0="WD-${instance}-"
solutionFile="outfiles/n060w8_0_6-2-9-9-0-8-1-3/OptSol-n060w8-62990813-"
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

java -jar validator.jar --sce datasets/n060w8/Sc-n060w8.txt --his datasets/n060w8/H0-n060w8-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n060w8_0_6-2-9-9-0-8-1-3/validatorOutput.txt 

exit 0;
