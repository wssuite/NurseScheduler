#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n050w4/Sc-n050w4.txt --his datasets/n050w4/H0-n050w4-0.txt --weeks datasets/n050w4/WD-n050w4-0.txt datasets/n050w4/WD-n050w4-4.txt datasets/n050w4/WD-n050w4-8.txt datasets/n050w4/WD-n050w4-7.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n050w4_0_0-4-8-7/seedTestsDemandingEvaluation/47-10-2-9 --rand 47 10 2 9 --timeout 113 --cus
