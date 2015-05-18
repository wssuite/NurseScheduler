#!/bin/bash

for entry in "outfiles/Competition/"*
do
	mv ${entry}/log.txt ${entry}/OptLog.txt 
	mv ${entry}/validatorOutput.txt ${entry}/OptValidatorOutput.txt
done
