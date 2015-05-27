#!/bin/bash

prefix="ess"
if [ ! -z "$1" ] then	prefix=$1
fi
echo $prefix

nbTests=30
stoOptions="stochastic_solver.txt"
geneOptios="generation_solver.txt"
evaOptions="evaluation_solver.txt"

#for entry in "sungrid"/*
#do
#	qsub ${entry} $1 $nbTests $stoOptions $geneOptios $evaOptions
#done
