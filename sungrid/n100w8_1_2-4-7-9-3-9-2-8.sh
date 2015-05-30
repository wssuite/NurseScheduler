#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n100w8 1 8 2 4 7 9 3 9 2 8 n100w8_1_2-4-7-9-3-9-2-8 $1 > outfiles/Competition/n100w8_1_2-4-7-9-3-9-2-8/${1}Log.txt

exit 0;
