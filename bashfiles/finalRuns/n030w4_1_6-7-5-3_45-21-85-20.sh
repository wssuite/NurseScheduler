#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n030w4/Sc-n030w4.txt --his datasets/n030w4/H0-n030w4-1.txt --weeks datasets/n030w4/WD-n030w4-6.txt datasets/n030w4/WD-n030w4-7.txt datasets/n030w4/WD-n030w4-5.txt datasets/n030w4/WD-n030w4-3.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n030w4_1_6-7-5-3/FinalRunAL/45-21-85-20/ --rand 45 21 85 20 --timeout 41 --cus
