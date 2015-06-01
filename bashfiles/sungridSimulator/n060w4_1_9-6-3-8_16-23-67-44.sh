#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n060w4/Sc-n060w4.txt --his datasets/n060w4/H0-n060w4-1.txt --weeks datasets/n060w4/WD-n060w4-9.txt datasets/n060w4/WD-n060w4-6.txt datasets/n060w4/WD-n060w4-3.txt datasets/n060w4/WD-n060w4-8.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n060w4_1_9-6-3-8/seedTestsDemandingEvaluation/16-23-67-44 --rand 16 23 67 44 --timeout 147 --cus
