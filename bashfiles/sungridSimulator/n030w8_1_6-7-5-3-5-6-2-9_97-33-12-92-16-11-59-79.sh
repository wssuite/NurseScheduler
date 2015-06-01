#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n030w8/Sc-n030w8.txt --his datasets/n030w8/H0-n030w8-1.txt --weeks datasets/n030w8/WD-n030w8-6.txt datasets/n030w8/WD-n030w8-7.txt datasets/n030w8/WD-n030w8-5.txt datasets/n030w8/WD-n030w8-3.txt datasets/n030w8/WD-n030w8-5.txt datasets/n030w8/WD-n030w8-6.txt datasets/n030w8/WD-n030w8-2.txt datasets/n030w8/WD-n030w8-9.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n030w8_1_6-7-5-3-5-6-2-9/seedTestsDemandingEvaluation/97-33-12-92-16-11-59-79 --rand 97 33 12 92 16 11 59 79 --timeout 45 --cus
