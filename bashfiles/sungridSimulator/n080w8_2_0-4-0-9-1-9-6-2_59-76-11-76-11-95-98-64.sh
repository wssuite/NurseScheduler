#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n080w8/Sc-n080w8.txt --his datasets/n080w8/H0-n080w8-2.txt --weeks datasets/n080w8/WD-n080w8-0.txt datasets/n080w8/WD-n080w8-4.txt datasets/n080w8/WD-n080w8-0.txt datasets/n080w8/WD-n080w8-9.txt datasets/n080w8/WD-n080w8-1.txt datasets/n080w8/WD-n080w8-9.txt datasets/n080w8/WD-n080w8-6.txt datasets/n080w8/WD-n080w8-2.txt --solver ./rosterDemandingEvaluation --runDir ./bin --outDir outfiles/n080w8_2_0-4-0-9-1-9-6-2/seedTestsDemandingEvaluation/59-76-11-76-11-95-98-64 --rand 59 76 11 76 11 95 98 64 --timeout 215 --cus
