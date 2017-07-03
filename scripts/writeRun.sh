# !/bin/sh

########################
# the first input is the instance name and the second one is the sequence of seeds
########################

# parse the input competition instance name
parse=(`echo $1 | tr '_' ' ' `)
arrayArgs=(`echo ${parse[*]} | tr '-' ' ' `)
instance=${arrayArgs[0]}
numHist=${arrayArgs[1]}
nbArgs=${#arrayArgs[@]}
weekList=${parse[2]}
nbWeeks=$((nbArgs - 2))

# get the parameters file from the second input
paramFile="paramfiles/${2}"
parseParam=(`echo $2 | tr '.' ' ' `)
paramName=${parseParam[0]}
seed=0

# display the inputs
echo "Instance: $1"
echo "History number: ${parse[1]}"
echo "Week list: ${parse[2]}"
echo "Parsed instance: ${arrayArgs[*]}"
echo "Parameter file: ${paramFile}"

# create the output directory if it does not exist
outputDir="outfiles/${paramName}/${1}"
if test ! -d "outfiles" ; then
	echo "Create output directory"
	mkdir "outfiles"
fi
if test ! -d "outfiles/${paramName}" ; then
	echo "Create deterministic output directory for parameters ${2}"
	mkdir "outfiles/${paramName}"
fi
if test ! -d "${outputDir}" ; then
	echo "Create this specific output directory"
	mkdir "${outputDir}"
fi

# create the files that are going to be used in the execution
scenarioFile="datasets/${instance}/Sc-${instance}.txt"
historyFile="datasets/${instance}/H0-${instance}-${numHist}.txt"

for ((i=2; i<$nbArgs; i++)); do
	demandFiles[$i-2]="datasets/${instance}/WD-${instance}-${arrayArgs[$i]}.txt"
done

for ((i=0; i<$nbWeeks; i++)); do
	solutionFiles[$i]="${outputDir}/sol-week${i}.txt"
done
validatorLog="${outputDir}/validator.txt"

# set the timeout depending on the operating system and number of nurses
if [[ "$OSTYPE" == "linux-gnu" ]]; then
	cpuMaxFor5Nurses=60
	cpuMaxPer10Nurses=60
elif [[ "$OSTYPE" == "darwin"* ]]; then
	cpuMaxFor5Nurses=60
	cpuMaxPer10Nurses=60
fi
nbTenNurses=${instance:1:2}
nbTenNurses=${nbTenNurses#0}
timeout=$((($nbWeeks)*($nbTenNurses)*$cpuMaxPer10Nurses+$cpuMaxFor5Nurses)) #43200 
if  [ "$2" == "optimality.txt" ] ; then
	timeout=86400
fi

# print the files that are going to be used during the execution
echo "Instance: ${instance}"
echo "Number of weeks: ${nbWeeks}"
echo "Scenario file: ${scenarioFile}"
echo "History file: ${historyFile}"
echo "Week files: ${demandFiles[*]}"
echo "Output directory: ${outputDir}"
echo "timeout = $timeout"
echo "seed = $seed"

# create the bashfile directory if it does not exist
if test ! -d "bashfiles" ; then
	echo "Create run directory"
	mkdir "bashfiles"
fi
if test ! -d "bashfiles/${paramName}" ; then
	echo "Create run directory for parameters ${paramName}"
	mkdir "bashfiles/${paramName}"
fi
bashfile="bashfiles/${paramName}/${1}.sh"

# write the instructions in the bashfile
echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q ${3}
#
# optimal script: run static scheduler" > ${bashfile}

echo "./bin/staticscheduler --dir datasets/ --instance ${instance} --weeks ${weekList} --his ${numHist} --param ${paramFile} --sol ${outputDir} --timeout ${timeout} 2> ${outputDir}/log.txt"  >> ${bashfile}
echo "java -jar validator.jar --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${solutionFiles[*]} > ${validatorLog}" >> ${bashfile}

chmod 755 "${bashfile}"
