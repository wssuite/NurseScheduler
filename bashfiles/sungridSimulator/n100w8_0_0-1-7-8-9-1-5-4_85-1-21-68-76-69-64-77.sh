#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n100w8/Sc-n100w8.txt --his datasets/n100w8/H0-n100w8-0.txt --weeks datasets/n100w8/WD-n100w8-0.txt datasets/n100w8/WD-n100w8-1.txt datasets/n100w8/WD-n100w8-7.txt datasets/n100w8/WD-n100w8-8.txt datasets/n100w8/WD-n100w8-9.txt datasets/n100w8/WD-n100w8-1.txt datasets/n100w8/WD-n100w8-5.txt datasets/n100w8/WD-n100w8-4.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n100w8_0_0-1-7-8-9-1-5-4/seedTestsDemandingEvaluation/85-1-21-68-76-69-64-77 --rand 85 1 21 68 76 69 64 77 --timeout 283 --cus
