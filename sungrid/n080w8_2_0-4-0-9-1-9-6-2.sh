#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n080w8 2 8 0 4 0 9 1 9 6 2 n080w8_2_0-4-0-9-1-9-6-2 $1 $2 $3 > outfiles/Competition/n080w8_2_0-4-0-9-1-9-6-2/${3}Log.txt

exit 0;
