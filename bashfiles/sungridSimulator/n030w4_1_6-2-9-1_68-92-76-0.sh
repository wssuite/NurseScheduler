#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n030w4/Sc-n030w4.txt --his datasets/n030w4/H0-n030w4-1.txt --weeks datasets/n030w4/WD-n030w4-6.txt datasets/n030w4/WD-n030w4-2.txt datasets/n030w4/WD-n030w4-9.txt datasets/n030w4/WD-n030w4-1.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n030w4_1_6-2-9-1/seedTestsDemandingEvaluation/68-92-76-0 --rand 68 92 76 0 --timeout 45 --cus
