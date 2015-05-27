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
./bin/optimalRoster n080w8 2 8 0 4 0 9 1 9 6 2 n080w8_2_0-4-0-9-1-9-6-2 $prefix $2 $3 $4 $5 > outfiles/Competition/n080w8_2_0-4-0-9-1-9-6-2/${1}Log.txt

instance=n080w8
weeksValue=(0 4 0 9 1 9 6 2 )

demand0="WD-${instance}-"
solutionFile="outfiles/Competition/n080w8_2_0-4-0-9-1-9-6-2/${1}sol-week"
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

java -jar validator.jar --sce datasets/n080w8/Sc-n080w8.txt --his datasets/n080w8/H0-n080w8-2.txt --weeks $weeks --sols $sols > outfiles/Competition/n080w8_2_0-4-0-9-1-9-6-2/${1}ValidatorOutput.txt 

exit 0;
