#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w4 1 4 9 6 3 8 n060w4_1_9-6-3-8 $1 $2 $3 > outfiles/Competition/n060w4_1_9-6-3-8/${3}Log.txt

exit 0;
