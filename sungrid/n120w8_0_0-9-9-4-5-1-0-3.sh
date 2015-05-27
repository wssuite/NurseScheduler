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
./bin/optimalRoster n120w8 0 8 0 9 9 4 5 1 0 3 n120w8_0_0-9-9-4-5-1-0-3 $prefix $2 $3 $4 $5 > outfiles/Competition/n120w8_0_0-9-9-4-5-1-0-3/${1}Log.txt

instance=n120w8
weeksValue=(0 9 9 4 5 1 0 3 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n120w8_0_0-9-9-4-5-1-0-3/${1}sol-week"
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

java -jar validator.jar --sce datasets/n120w8/Sc-n120w8.txt --his datasets/n120w8/H0-n120w8-0.txt --weeks $weeks --sols $sols > outfiles/Competition/n120w8_0_0-9-9-4-5-1-0-3/${1}ValidatorOutput.txt 

exit 0;
