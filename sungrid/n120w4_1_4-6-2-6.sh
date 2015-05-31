#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n120w4 1 4 4 6 2 6 n120w4_1_4-6-2-6 $1 > outfiles/Competition/n120w4_1_4-6-2-6/${1}Log.txt

exit 0;
