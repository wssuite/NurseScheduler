#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n120w4/Sc-n120w4.txt --his datasets/n120w4/H0-n120w4-1.txt --weeks datasets/n120w4/WD-n120w4-5.txt datasets/n120w4/WD-n120w4-6.txt datasets/n120w4/WD-n120w4-9.txt datasets/n120w4/WD-n120w4-8.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n120w4_1_5-6-9-8/FinalRun/99-77-89-18/ --rand 99 77 89 18 --timeout 320 --cus
