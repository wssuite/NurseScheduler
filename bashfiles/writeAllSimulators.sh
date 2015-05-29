# !/bin/sh

# parse the input competition instance name
CURDIR=`pwd`
cd outfiles/Competition
for DIR in * 
do 
	cd ${CURDIR}
	bashfiles/writeSimulatorWithRandomSeeds.sh ${DIR}
done