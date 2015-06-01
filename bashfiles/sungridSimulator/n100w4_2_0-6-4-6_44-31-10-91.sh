#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n100w4/Sc-n100w4.txt --his datasets/n100w4/H0-n100w4-2.txt --weeks datasets/n100w4/WD-n100w4-0.txt datasets/n100w4/WD-n100w4-6.txt datasets/n100w4/WD-n100w4-4.txt datasets/n100w4/WD-n100w4-6.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n100w4_2_0-6-4-6/seedTestsDemandingEvaluation/44-31-10-91 --rand 44 31 10 91 --timeout 283 --cus
