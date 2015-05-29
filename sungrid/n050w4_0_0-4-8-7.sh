#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n050w4 0 4 0 4 8 7 n050w4_0_0-4-8-7 $1 $2 $3 > outfiles/Competition/n050w4_0_0-4-8-7/${3}Log.txt

exit 0;
