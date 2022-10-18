# !/bin/sh

# create the scripts for all the deterministic runs at once
CURDIR=`pwd`
INSTANCES=()

for inputFile in datasets/INRC/*.txt
do
  echo ${inputFile}
  INSTANCES+=(${inputFile})
done

mkdir -p "bashfiles"
mkdir -p "bashfiles/INRC"
mkdir -p "outfiles"
mkdir -p "outfiles/INRC"


for INST in ${INSTANCES[*]}
do
  filename="${INST##*/}"
  dirname="${INST%$filename}"
  instname="${filename%.*}"
  bashfile="bashfiles/INRC/${instname}.sh"
  validatorFile="outfiles/INRC/${instname}_validator.log"
  logFile="outfiles/INRC/${instname}.log"
  echo $bashfile

  echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#
# optimal script: run static scheduler on INRC datasets" > ${bashfile}
	echo "./bin/inrc ${dirname} ${filename} > ${logFile}"  >> ${bashfile}
  echo "java -jar validatorINRC.jar datasets/INRC/instances_xml/${instname}.xml outfiles/INRC/${instname}.xml > ${validatorFile}" >> ${bashfile}
  chmod 755 "${bashfile}"
done
