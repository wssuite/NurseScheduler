#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch optimal solver and then the validator

./bin/optimalRoster n120w8 0 8 0 9 9 4 5 1 0 3 n120w8_0_0-9-9-4-5-1-0-3 $1 $2 $3 > outfiles/Competition/n120w8_0_0-9-9-4-5-1-0-3/${3}Log.txt

exit 0;
