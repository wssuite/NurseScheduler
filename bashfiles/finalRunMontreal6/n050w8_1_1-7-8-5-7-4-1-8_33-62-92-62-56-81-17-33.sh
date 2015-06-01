#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n050w8/Sc-n050w8.txt --his datasets/n050w8/H0-n050w8-1.txt --weeks datasets/n050w8/WD-n050w8-1.txt datasets/n050w8/WD-n050w8-7.txt datasets/n050w8/WD-n050w8-8.txt datasets/n050w8/WD-n050w8-5.txt datasets/n050w8/WD-n050w8-7.txt datasets/n050w8/WD-n050w8-4.txt datasets/n050w8/WD-n050w8-1.txt datasets/n050w8/WD-n050w8-8.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n050w8_1_1-7-8-5-7-4-1-8/FinalRun/33-62-92-62-56-81-17-33/ --rand 33 62 92 62 56 81 17 33 --timeout 103 --cus
