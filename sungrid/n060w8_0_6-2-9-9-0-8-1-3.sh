#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n060w8 0 8 6 2 9 9 0 8 1 3 n060w8_0_6-2-9-9-0-8-1-3 $1 $2 $3 > outfiles/Competition/n060w8_0_6-2-9-9-0-8-1-3/${3}Log.txt

exit 0;
