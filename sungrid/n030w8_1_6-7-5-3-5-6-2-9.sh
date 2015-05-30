#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n030w8 1 8 6 7 5 3 5 6 2 9 n030w8_1_6-7-5-3-5-6-2-9 $1 > outfiles/Competition/n030w8_1_6-7-5-3-5-6-2-9/${1}Log.txt

exit 0;
