#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
java -jar Simulator.jar  --sce datasets/n100w4/Sc-n100w4.txt --his datasets/n100w4/H0-n100w4-0.txt --weeks datasets/n100w4/WD-n100w4-1.txt datasets/n100w4/WD-n100w4-1.txt datasets/n100w4/WD-n100w4-0.txt datasets/n100w4/WD-n100w4-8.txt --solver ./bin/roster --runDir ./ --outDir outfiles/n100w4_0_1-1-0-8/
