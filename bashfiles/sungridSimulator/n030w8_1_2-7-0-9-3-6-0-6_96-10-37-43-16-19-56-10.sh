#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n030w8/Sc-n030w8.txt --his datasets/n030w8/H0-n030w8-1.txt --weeks datasets/n030w8/WD-n030w8-2.txt datasets/n030w8/WD-n030w8-7.txt datasets/n030w8/WD-n030w8-0.txt datasets/n030w8/WD-n030w8-9.txt datasets/n030w8/WD-n030w8-3.txt datasets/n030w8/WD-n030w8-6.txt datasets/n030w8/WD-n030w8-0.txt datasets/n030w8/WD-n030w8-6.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n030w8_1_2-7-0-9-3-6-0-6/seedTestsDemandingEvaluation/96-10-37-43-16-19-56-10 --rand 96 10 37 43 16 19 56 10 --timeout 45 --cus
