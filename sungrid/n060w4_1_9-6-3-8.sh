#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n060w4 1 4 9 6 3 8 n060w4_1_9-6-3-8 $1 > outfiles/Competition/n060w4_1_9-6-3-8/${1}Log.txt

exit 0;
