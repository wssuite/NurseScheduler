#!/bin/bash

for entry in "sungrid"/*
do
#seeds
	for i in {1..30}
	do
#nb evaluation
		for n in {1..10}
		do
			qsub ${entry} $i $n $1
		done
	done
done
