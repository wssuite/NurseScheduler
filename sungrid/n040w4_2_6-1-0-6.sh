#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n040w4 2 4 6 1 0 6 n040w4_2_6-1-0-6 $1 > outfiles/Competition/n040w4_2_6-1-0-6/${1}Log.txt

exit 0;
