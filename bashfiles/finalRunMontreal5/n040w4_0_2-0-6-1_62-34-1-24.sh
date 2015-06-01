#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n040w4/Sc-n040w4.txt --his datasets/n040w4/H0-n040w4-0.txt --weeks datasets/n040w4/WD-n040w4-2.txt datasets/n040w4/WD-n040w4-0.txt datasets/n040w4/WD-n040w4-6.txt datasets/n040w4/WD-n040w4-1.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n040w4_0_2-0-6-1/FinalRun/62-34-1-24/ --rand 62 34 1 24 --timeout 72 --cus
