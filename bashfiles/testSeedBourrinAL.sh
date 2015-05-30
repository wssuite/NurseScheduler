# !/bin/sh

# parse the input competition instance name
CURDIR=`pwd`
NB_TESTS_PER_INSTANCE=1

BASHDIR="bashfiles/sungridSimulator/${RANDOM}"
if test ! -d "bashfiles/sungridSimulator" ; then
	echo "Create output directory bashfiles/sungridSimulator"
	mkdir "bashfiles/sungridSimulator"
fi
if test ! -d "${outputDir}" ; then
	echo "Create output directory ${BASHDIR}"
	mkdir "${BASHDIR}"
fi

cd outfiles/Competition
for DIR in * 
do 
	cd ${CURDIR}
	for ((i=0; i<$NB_TESTS_PER_INSTANCE; i++)); do
		bashfiles/writeSimulatorWithRandomSeedsAL.sh ${DIR} ${BASHDIR}
	done
done

for BASHFILE in ${BASHDIR}/*
do
  qsub ${BASHFILE}
done

