#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n050w8/Sc-n050w8.txt --his datasets/n050w8/H0-n050w8-1.txt --weeks datasets/n050w8/WD-n050w8-9.txt datasets/n050w8/WD-n050w8-7.txt datasets/n050w8/WD-n050w8-5.txt datasets/n050w8/WD-n050w8-3.txt datasets/n050w8/WD-n050w8-8.txt datasets/n050w8/WD-n050w8-8.txt datasets/n050w8/WD-n050w8-3.txt datasets/n050w8/WD-n050w8-1.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n050w8_1_9-7-5-3-8-8-3-1/seedTestsDemandingEvaluation/1-25-36-56-31-29-64-4 --rand 1 25 36 56 31 29 64 4 --timeout 113 --cus
