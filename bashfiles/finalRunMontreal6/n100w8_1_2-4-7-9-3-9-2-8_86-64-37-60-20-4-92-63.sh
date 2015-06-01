#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n100w8/Sc-n100w8.txt --his datasets/n100w8/H0-n100w8-1.txt --weeks datasets/n100w8/WD-n100w8-2.txt datasets/n100w8/WD-n100w8-4.txt datasets/n100w8/WD-n100w8-7.txt datasets/n100w8/WD-n100w8-9.txt datasets/n100w8/WD-n100w8-3.txt datasets/n100w8/WD-n100w8-9.txt datasets/n100w8/WD-n100w8-2.txt datasets/n100w8/WD-n100w8-8.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n100w8_1_2-4-7-9-3-9-2-8/FinalRun/86-64-37-60-20-4-92-63/ --rand 86 64 37 60 20 4 92 63 --timeout 258 --cus
