#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n100w8 0 8 0 1 7 8 9 1 5 4 n100w8_0_0-1-7-8-9-1-5-4 $1 > outfiles/Competition/n100w8_0_0-1-7-8-9-1-5-4/${1}Log.txt

exit 0;
