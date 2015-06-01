#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n060w8/Sc-n060w8.txt --his datasets/n060w8/H0-n060w8-2.txt --weeks datasets/n060w8/WD-n060w8-1.txt datasets/n060w8/WD-n060w8-0.txt datasets/n060w8/WD-n060w8-3.txt datasets/n060w8/WD-n060w8-4.txt datasets/n060w8/WD-n060w8-0.txt datasets/n060w8/WD-n060w8-3.txt datasets/n060w8/WD-n060w8-9.txt datasets/n060w8/WD-n060w8-1.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n060w8_2_1-0-3-4-0-3-9-1/seedTestsDemandingEvaluation/43-58-36-91-83-2-2-8 --rand 43 58 36 91 83 2 2 8 --timeout 147 --cus
