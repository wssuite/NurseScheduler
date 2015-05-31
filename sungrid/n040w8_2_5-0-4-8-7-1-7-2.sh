#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n040w8 2 8 5 0 4 8 7 1 7 2 n040w8_2_5-0-4-8-7-1-7-2 $1 > outfiles/Competition/n040w8_2_5-0-4-8-7-1-7-2/${1}Log.txt

exit 0;
