# !/bin/sh

# parse the input competition instance name
CURDIR=`pwd`
NB_TESTS_PER_INSTANCE=50

for ((j=0; j<$NB_TESTS_PER_INSTANCE; j++)); do
BASHDIR="bashfiles/sungridSimulator/Flexible${RANDOM}"
if test ! -d "bashfiles/sungridSimulator" ; then
	echo "Create output directory bashfiles/sungridSimulator"
	mkdir "bashfiles/sungridSimulator"
fi
if test ! -d "${BASHDIR}" ; then
	echo "Create output directory ${BASHDIR}"
	mkdir "${BASHDIR}"
fi

cd outfiles/Competition
for DIR in * 
do 
	cd ${CURDIR}
	sh bashfiles/writeSimulatorWithRandomSeedsAL.sh ${DIR} ${BASHDIR}
done

for BASHFILE in ${BASHDIR}/*
do
  qsub ${BASHFILE}
done

done
