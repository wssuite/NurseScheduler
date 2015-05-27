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
./bin/optimalRoster n100w4 0 4 1 1 0 8 n100w4_0_1-1-0-8 $prefix $2 $3 $4 $5 > outfiles/Competition/n100w4_0_1-1-0-8/${1}Log.txt

instance=n100w4
weeksValue=(1 1 0 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n100w4_0_1-1-0-8/${1}sol-week"
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

java -jar validator.jar --sce datasets/n100w4/Sc-n100w4.txt --his datasets/n100w4/H0-n100w4-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n100w4_0_1-1-0-8/${1}ValidatorOutput.txt 

exit 0;
