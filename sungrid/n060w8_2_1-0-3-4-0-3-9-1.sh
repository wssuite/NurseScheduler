#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w8 2 8 1 0 3 4 0 3 9 1 n060w8_2_1-0-3-4-0-3-9-1 $1 $2 $3 > outfiles/Competition/n060w8_2_1-0-3-4-0-3-9-1/${3}Log.txt

exit 0;
