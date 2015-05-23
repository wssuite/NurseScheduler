#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n080w4/Sc-n080w4.txt --his datasets/n080w4/H0-n080w4-2.txt --weeks datasets/n080w4/WD-n080w4-4.txt datasets/n080w4/WD-n080w4-3.txt datasets/n080w4/WD-n080w4-3.txt datasets/n080w4/WD-n080w4-3.txt --solver ./roster --runDir ./bin --outDir outfiles/n080w4_2_4-3-3-3/WithNoEvaluation/
