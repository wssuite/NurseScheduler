#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n100w4 0 4 1 1 0 8 n100w4_0_1-1-0-8 $1 $2 $3 > outfiles/Competition/n100w4_0_1-1-0-8/${3}Log.txt

exit 0;
