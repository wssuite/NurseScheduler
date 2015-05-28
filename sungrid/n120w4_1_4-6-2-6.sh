#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n120w4 1 4 4 6 2 6 n120w4_1_4-6-2-6 $1 $2 $3 > outfiles/Competition/n120w4_1_4-6-2-6/${3}Log.txt

exit 0;
