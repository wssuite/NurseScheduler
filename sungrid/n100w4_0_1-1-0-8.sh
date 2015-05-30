#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q osiris
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n100w4 0 4 1 1 0 8 n100w4_0_1-1-0-8 $1 > outfiles/Competition/n100w4_0_1-1-0-8/${1}Log.txt

exit 0;
