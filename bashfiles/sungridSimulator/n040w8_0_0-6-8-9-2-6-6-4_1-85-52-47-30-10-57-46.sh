#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n040w8/Sc-n040w8.txt --his datasets/n040w8/H0-n040w8-0.txt --weeks datasets/n040w8/WD-n040w8-0.txt datasets/n040w8/WD-n040w8-6.txt datasets/n040w8/WD-n040w8-8.txt datasets/n040w8/WD-n040w8-9.txt datasets/n040w8/WD-n040w8-2.txt datasets/n040w8/WD-n040w8-6.txt datasets/n040w8/WD-n040w8-6.txt datasets/n040w8/WD-n040w8-4.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n040w8_0_0-6-8-9-2-6-6-4/seedTestsDemandingEvaluation/1-85-52-47-30-10-57-46 --rand 1 85 52 47 30 10 57 46 --timeout 79 --cus
