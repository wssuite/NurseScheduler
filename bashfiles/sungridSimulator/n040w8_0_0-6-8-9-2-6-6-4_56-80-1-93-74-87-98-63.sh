#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n040w8/Sc-n040w8.txt --his datasets/n040w8/H0-n040w8-0.txt --weeks datasets/n040w8/WD-n040w8-0.txt datasets/n040w8/WD-n040w8-6.txt datasets/n040w8/WD-n040w8-8.txt datasets/n040w8/WD-n040w8-9.txt datasets/n040w8/WD-n040w8-2.txt datasets/n040w8/WD-n040w8-6.txt datasets/n040w8/WD-n040w8-6.txt datasets/n040w8/WD-n040w8-4.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n040w8_0_0-6-8-9-2-6-6-4/seedTestsDemandingEvaluation/56-80-1-93-74-87-98-63 --rand 56 80 1 93 74 87 98 63 --timeout 79 --cus
