#!/bin/bash

nbTests=30
#stoOptions="stochastic_solver.txt"
geneOptios="generation_solver.txt"
evaOptions="evaluation_solver.txt"

for entry in "sungrid"/*
do
for n in {1..10}
do
	stoOptions="stochastic_solver$n.txt"
	qsub ${entry} $1 $nbTests $stoOptions $geneOptios $evaOptions
done
done
