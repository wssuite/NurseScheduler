#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n060w4 1 4 6 1 1 5 n060w4_1_6-1-1-5 $1 > outfiles/Competition/n060w4_1_6-1-1-5/${1}Log.txt

exit 0;
