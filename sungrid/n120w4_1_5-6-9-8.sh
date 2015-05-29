#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n120w4 1 4 5 6 9 8 n120w4_1_5-6-9-8 $1 $2 $3 > outfiles/Competition/n120w4_1_5-6-9-8/${3}Log.txt

exit 0;
