#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

prefix="-1"
if [ ! -z "$1" ]
then	prefix=$1
fi
./bin/optimalRoster n060w4 1 4 9 6 3 8 n060w4_1_9-6-3-8 $prefix $2 $3 $4 $5 > outfiles/Competition/n060w4_1_9-6-3-8/${1}Log.txt

instance=n060w4
weeksValue=(9 6 3 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n060w4_1_9-6-3-8/${1}sol-week"
weeks=""
sols=""
i=0

for var in ${weeksValue[*]}
do
demand[$i]="datasets/${instance}/${demand0}${var}.txt"
weeks="${weeks} ${demand[$i]}"
solution[$i]="${solutionFile}${i}.txt"
sols="${sols} ${solution[$i]}"
((i++))
done

java -jar validator.jar --sce datasets/n060w4/Sc-n060w4.txt --his datasets/n060w4/H0-n060w4-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n060w4_1_9-6-3-8/${1}ValidatorOutput.txt 

exit 0;
