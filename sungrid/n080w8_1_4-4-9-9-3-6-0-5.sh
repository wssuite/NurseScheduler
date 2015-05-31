#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/rankingRoster n080w8 1 8 4 4 9 9 3 6 0 5 n080w8_1_4-4-9-9-3-6-0-5 $1 > outfiles/Competition/n080w8_1_4-4-9-9-3-6-0-5/${1}Log.txt

exit 0;
