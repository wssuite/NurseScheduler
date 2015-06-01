#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n060w8/Sc-n060w8.txt --his datasets/n060w8/H0-n060w8-0.txt --weeks datasets/n060w8/WD-n060w8-6.txt datasets/n060w8/WD-n060w8-2.txt datasets/n060w8/WD-n060w8-9.txt datasets/n060w8/WD-n060w8-9.txt datasets/n060w8/WD-n060w8-0.txt datasets/n060w8/WD-n060w8-8.txt datasets/n060w8/WD-n060w8-1.txt datasets/n060w8/WD-n060w8-3.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/Competition/n060w8_0_6-2-9-9-0-8-1-3/FinalRun/96-24-68-88-53-91-65-36/ --rand 96 24 68 88 53 91 65 36 --timeout 134 --cus
