#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n080w4 2 4 6 0 4 8 n080w4_2_6-0-4-8 $1 > outfiles/Competition/n080w4_2_6-0-4-8/${1}Log.txt

exit 0;
