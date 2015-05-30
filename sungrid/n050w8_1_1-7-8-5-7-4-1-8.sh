#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n050w8 1 8 1 7 8 5 7 4 1 8 n050w8_1_1-7-8-5-7-4-1-8 $1 > outfiles/Competition/n050w8_1_1-7-8-5-7-4-1-8/${1}Log.txt

exit 0;
