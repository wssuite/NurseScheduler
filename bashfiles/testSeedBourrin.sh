# !/bin/sh

# parse the input competition instance name
CURDIR=`pwd`
NB_TESTS_PER_INSTANCE=10

rm bashfiles/sungridSimulator/*
cd outfiles/Competition
for DIR in * 
do 
	cd ${CURDIR}
	for ((i=0; i<$NB_TESTS_PER_INSTANCE; i++)); do
		bashfiles/writeSimulatorWithRandomSeeds.sh ${DIR}
	done
done

for BASHFILE in "bashfiles/sungridSimulator"/*
do
  qsub ${BASHFILE}
done

