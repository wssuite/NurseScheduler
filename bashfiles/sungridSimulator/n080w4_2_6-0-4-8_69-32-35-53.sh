#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n080w4/Sc-n080w4.txt --his datasets/n080w4/H0-n080w4-2.txt --weeks datasets/n080w4/WD-n080w4-6.txt datasets/n080w4/WD-n080w4-0.txt datasets/n080w4/WD-n080w4-4.txt datasets/n080w4/WD-n080w4-8.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n080w4_2_6-0-4-8/seedTestsDemandingEvaluation/69-32-35-53 --rand 69 32 35 53 --timeout 215 --cus
