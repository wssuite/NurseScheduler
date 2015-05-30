#!/bin/bash

for entry in "sungrid"/*
do
	for n in {1..30}
	do
		qsub ${entry} $1
	done
done
