#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n120w8 1 8 7 2 6 4 5 2 0 2 n120w8_1_7-2-6-4-5-2-0-2 $1 > outfiles/Competition/n120w8_1_7-2-6-4-5-2-0-2/${1}Log.txt

exit 0;
