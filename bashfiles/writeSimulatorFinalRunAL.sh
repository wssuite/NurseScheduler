# !/bin/sh

########################
# the first input is the instance name and the second one is the sequence of seeds
########################

# parse the input competition instance name
echo "Competition instance: $1"
parse=(`echo $1 | tr '_' ' ' `)
arrayArgs=(`echo ${parse[*]} | tr '-' ' ' `)
echo "${arrayArgs[*]}"
instance=${arrayArgs[0]}
numHist=${arrayArgs[1]}
nbArgs=${#arrayArgs[@]}
nbWeeks=$((nbArgs - 2))

# get the seeds form the second input
seeds=(`echo $2 | tr '-' ' ' `)

# create the output directory if it does not exist
outputDir="outfiles/Competition/${1}/"
if test ! -d "outfiles" ; then
	echo "Create output directory"
	mkdir "outfiles"
fi
if test ! -d "${outputDir}" ; then
	echo "Create output directory"
	mkdir "${outputDir}"
fi
outputDir="outfiles/Competition/${1}/FinalRun/"
if test ! -d "${outputDir}" ; then
	echo "Create output directory"
	mkdir "${outputDir}"
fi

# create the files that are going to be used in the execution
scenarioFile="datasets/${instance}/Sc-${instance}.txt"
historyFile="datasets/${instance}/H0-${instance}-${numHist}.txt"

for ((i=2; i<$nbArgs; i++)); do
	demandFiles[$i-2]="datasets/${instance}/WD-${instance}-${arrayArgs[$i]}.txt"
done

# set the timeout depending on the operating system and number of nurses
nbNurses=${instance:1:3}
nbNurses=${nbNurses#0}
timeout=$((($nbNurses-20)*3+10))

# print the files that are going to be used during the execution
echo "Instance: ${instance}"
echo "Number of weeks: ${nbWeeks}"
echo "Scenario file: ${scenarioFile}"
echo "History file: ${historyFile}"
echo "Week files: ${demandFiles[*]}"
echo "Output directory: ${outputDir}"
echo "timeout = $timeout"
echo "seeds = ${seeds[*]}"

if test ! -d "bashfiles/finalRuns" ; then
	echo "Create run directory"
	mkdir "bashfiles/finalRuns"
fi
bashfile="bashfiles/finalRuns/${1}.sh"

echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator" > ${bashfile}

echo "java -jar Simulator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --solver ./roster --runDir ./bin --outDir ${outputDir} --rand ${seeds[*]} --timeout ${timeout} --cus"  >> ${bashfile}

chmod 755 "${bashfile}"
#echo "java -jar validator.jar --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${solutionFiles[*]} > ${validatorLog}" >> ${sungridfile}
