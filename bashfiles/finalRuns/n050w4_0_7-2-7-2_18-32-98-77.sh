#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n050w4/Sc-n050w4.txt --his datasets/n050w4/H0-n050w4-0.txt --weeks datasets/n050w4/WD-n050w4-7.txt datasets/n050w4/WD-n050w4-2.txt datasets/n050w4/WD-n050w4-7.txt datasets/n050w4/WD-n050w4-2.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n050w4_0_7-2-7-2/FinalRunAL/18-32-98-77/ --rand 18 32 98 77 --timeout 103 --cus
