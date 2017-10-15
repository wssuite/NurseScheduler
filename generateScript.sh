# !/bin/sh

# $1: contains the description of the instance: n005w4_1-2-3-3_0 (data_weeks_history)
# $2: contains the seeds: 22-36-96-5. If not defined, generate random seeds

# parse the input competition instance name
echo "Competition instance: $1"
parse=(`echo $1 | tr '_' ' ' `)
arrayArgs=(`echo ${parse[*]} | tr '-' ' ' `)
echo "${arrayArgs[*]}"
instance=${arrayArgs[0]}
numHist=${arrayArgs[1]}
nbArgs=${#arrayArgs[@]}
nbWeeks=$((nbArgs - 2))

# create the files that are going to be used in the execution
scenarioFile="datasets/${instance}/Sc-${instance}.txt"
historyFile="datasets/${instance}/H0-${instance}-${numHist}.txt"

for ((i=2; i<$nbArgs; i++)); do
	demandFiles[$i-2]="datasets/${instance}/WD-${instance}-${arrayArgs[$i]}.txt"
done

allseeds=$2
if [ -z $2 ]
then
	# create the random seeds
	allseeds=$(($RANDOM%100))
	for ((i=1; i<$nbWeeks; i++))
	do
		allseeds+="-$(($RANDOM%100))"
	done
fi
# get the seeds
seeds=(`echo $allseeds | tr '-' ' ' `)

# create the root of output directory if it does not exist
outputDir="outfiles/${1}/$allseeds"
echo "Create output directory: ${outputDir}"
mkdir -p "${outputDir}"

# set the solutions vector
sols="${outputDir}/sol-week0.txt"
for ((i=1; i<$nbWeeks; i++)); do
	sols+=" ${outputDir}/sol-week${i}.txt"
done

# set the timeout depending on the operating system and number of nurses
nbNurses=${instance:1:3}
nbNurses=${nbNurses#0}
timeout=$((($nbNurses-20)*3+10))
if [ $timeout -lt 5 ]
then
	timeout=5
fi

# print the files that are going to be used during the execution
echo "Instance: ${instance}"
echo "Number of weeks: ${nbWeeks}"
echo "Scenario file: ${scenarioFile}"
echo "History file: ${historyFile}"
echo "Week files: ${demandFiles[*]}"
echo "Output directory: ${outputDir}"
echo "timeout = $timeout"
echo "seeds = ${seeds[*]}"
echo "solutions = ${sols}"
echo "log validator = ${outputDir}/validator.txt"

# script file
scriptfile="${1}_${allseeds}.sh"
echo "scriptfile=$scriptfile"

echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
" > ${scriptfile}

echo "java -jar Simulator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --solver ./roster --runDir ./bin --outDir ${outputDir} --rand ${seeds[*]} --timeout ${timeout} --cus"  >> ${scriptfile}

echo "java -jar validator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${sols} > ${outputDir}/validator.txt"  >> ${scriptfile}

chmod 755 "${scriptfile}"
