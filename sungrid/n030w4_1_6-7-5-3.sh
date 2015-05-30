#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n030w4 1 4 6 7 5 3 n030w4_1_6-7-5-3 $1 > outfiles/Competition/n030w4_1_6-7-5-3/${1}Log.txt

exit 0;
