#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n080w8/Sc-n080w8.txt --his datasets/n080w8/H0-n080w8-1.txt --weeks datasets/n080w8/WD-n080w8-4.txt datasets/n080w8/WD-n080w8-4.txt datasets/n080w8/WD-n080w8-9.txt datasets/n080w8/WD-n080w8-9.txt datasets/n080w8/WD-n080w8-3.txt datasets/n080w8/WD-n080w8-6.txt datasets/n080w8/WD-n080w8-0.txt datasets/n080w8/WD-n080w8-5.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n080w8_1_4-4-9-9-3-6-0-5/seedTestsDemandingEvaluation/93-91-1-97-80-17-35-52 --rand 93 91 1 97 80 17 35 52 --timeout 215 --cus
