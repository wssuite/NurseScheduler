#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n040w8 0 8 0 6 8 9 2 6 6 4 n040w8_0_0-6-8-9-2-6-6-4 $1 $2 $3 > outfiles/Competition/n040w8_0_0-6-8-9-2-6-6-4/${3}Log.txt

exit 0;
