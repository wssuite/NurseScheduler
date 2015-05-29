#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n050w4 0 4 7 2 7 2 n050w4_0_7-2-7-2 $1 $2 $3 > outfiles/Competition/n050w4_0_7-2-7-2/${3}Log.txt

exit 0;
