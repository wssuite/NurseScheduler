#!/bin/bash

entry="sungrid/n120w8_1_7-2-6-4-5-2-0-2.sh"
#for entry in "sungrid"/*
#do
#seeds
	for i in {1..30}
	do
#nb evaluation
		for n in {1..10}
		do
			qsub ${entry} $i $n $1
		done
	done
#done
