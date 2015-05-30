#!/bin/bash

#entry="sungrid/n030w8_1_2-7-0-9-3-6-0-6.sh"
#entry="sungrid/n060w4_1_9-6-3-8.sh"
entry="sungrid/n100w4_0_1-1-0-8.sh"
#entry="sungrid/n080w8_1_4-4-9-9-3-6-0-5.sh"

#./bin/nostratMeanRoster n080w8 1 8 4 4 9 9 3 6 0 5 n080w8_1_4-4-9-9-3-6-0-5 $s 1 nostrat_mean_ 
#./bin/nostratMeanRoster n100w4 0 4 1 1 0 8 n100w4_0_1-1-0-8 $s 1 nostrat_mean_
#./bin/nostratMeanRoster n060w4 1 4 9 6 3 8 n060w4_1_9-6-3-8 $s 1 nostrat_ mean_
#./bin/nostratMeanRoster n030w8 1 8 2 7 0 9 3 6 0 6 n030w8_1_2-7-0-9-3-6-0-6 $s 1 nostrat_mean_

#for entry in "sungrid"/*
#do
#seeds
	for i in {1..30}
	do
#nb evaluation
		for n in {0..10}
		do
			sh ${entry} $i $n $1
		done
	done
#done
