#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n100w4 2 4 0 6 4 6 n100w4_2_0-6-4-6 $1 $2 $3 > outfiles/Competition/n100w4_2_0-6-4-6/${3}Log.txt

exit 0;
