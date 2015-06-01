#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n040w4/Sc-n040w4.txt --his datasets/n040w4/H0-n040w4-0.txt --weeks datasets/n040w4/WD-n040w4-2.txt datasets/n040w4/WD-n040w4-0.txt datasets/n040w4/WD-n040w4-6.txt datasets/n040w4/WD-n040w4-1.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n040w4_0_2-0-6-1/seedTestsDemandingEvaluation/0-87-67-46 --rand 0 87 67 46 --timeout 79 --cus
