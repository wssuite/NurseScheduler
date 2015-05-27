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
./bin/optimalRoster n080w4 2 4 6 0 4 8 n080w4_2_6-0-4-8 $prefix $2 $3 $4 $5 > outfiles/Competition/n080w4_2_6-0-4-8/${1}Log.txt

instance=n080w4
weeksValue=(6 0 4 8 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n080w4_2_6-0-4-8/${1}sol-week"
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

java -jar validator.jar --sce datasets/n080w4/Sc-n080w4.txt --his datasets/n080w4/H0-n080w4-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n080w4_2_6-0-4-8/${1}ValidatorOutput.txt 

exit 0;
