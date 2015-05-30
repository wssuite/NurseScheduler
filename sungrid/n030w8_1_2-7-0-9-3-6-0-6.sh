#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n030w8 1 8 2 7 0 9 3 6 0 6 n030w8_1_2-7-0-9-3-6-0-6 $1 > outfiles/Competition/n030w8_1_2-7-0-9-3-6-0-6/${1}Log.txt

exit 0;
