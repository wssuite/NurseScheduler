#!/bin/bash
prefix="Off"
for entry in "outfiles/Competition/"*
do
	cd ${entry}
	mv log.txt ${prefix}Log.txt 
	mv validatorOutput.txt ${prefix}ValidatorOutput.txt
	for entry2 in "OptSol"*
	do
		mv ${entry2} ${prefix}${entry2}
	done
	cd ../../..
done
