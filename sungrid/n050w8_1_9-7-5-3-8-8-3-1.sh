#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n050w8 1 8 9 7 5 3 8 8 3 1 n050w8_1_9-7-5-3-8-8-3-1 $1 > outfiles/Competition/n050w8_1_9-7-5-3-8-8-3-1/${1}Log.txt

exit 0;
