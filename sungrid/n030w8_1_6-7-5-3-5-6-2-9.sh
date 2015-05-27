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
./bin/optimalRoster n030w8 1 8 6 7 5 3 5 6 2 9 n030w8_1_6-7-5-3-5-6-2-9 $prefix $2 $3 $4 $5 > outfiles/Competition/n030w8_1_6-7-5-3-5-6-2-9/${1}Log.txt

instance=n030w8
weeksValue=(6 7 5 3 5 6 2 9 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n030w8_1_6-7-5-3-5-6-2-9/${1}sol-week"
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

java -jar validator.jar --sce datasets/n030w8/Sc-n030w8.txt --his datasets/n030w8/H0-n030w8-1.txt --weeks $weeks --sols $sols > outfiles/Competition/n030w8_1_6-7-5-3-5-6-2-9/${1}ValidatorOutput.txt 

exit 0;
