#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n030w8/Sc-n030w8.txt --his datasets/n030w8/H0-n030w8-1.txt --weeks datasets/n030w8/WD-n030w8-6.txt datasets/n030w8/WD-n030w8-7.txt datasets/n030w8/WD-n030w8-5.txt datasets/n030w8/WD-n030w8-3.txt datasets/n030w8/WD-n030w8-5.txt datasets/n030w8/WD-n030w8-6.txt datasets/n030w8/WD-n030w8-2.txt datasets/n030w8/WD-n030w8-9.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n030w8_1_6-7-5-3-5-6-2-9/FinalRun/73-57-84-45-83-2-62-87/ --rand 73 57 84 45 83 2 62 87 --timeout 41 --cus
