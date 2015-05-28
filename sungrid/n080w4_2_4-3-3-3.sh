#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n080w4 2 4 4 3 3 3 n080w4_2_4-3-3-3 $1 $2 $3 > outfiles/Competition/n080w4_2_4-3-3-3/${3}Log.txt

exit 0;
