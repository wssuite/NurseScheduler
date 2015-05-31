#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n040w4 0 4 2 0 6 1 n040w4_0_2-0-6-1 $1 > outfiles/Competition/n040w4_0_2-0-6-1/${1}Log.txt

exit 0;
