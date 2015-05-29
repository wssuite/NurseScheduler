#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n040w4 2 4 6 1 0 6 n040w4_2_6-1-0-6 $1 $2 $3 > outfiles/Competition/n040w4_2_6-1-0-6/${3}Log.txt

exit 0;
