# !/bin/sh

########################
# the first input is the instance name and the second one is the sequence of seeds
########################

function printBashUsage {
  echo "This script will run the simulator and then the validator."
  echo "Usage:"
  echo "-h | --help: display this message"
  echo "-i | --instance: instance to test (must follow the pattern (data_weeks_history)). Default: n030w4_1-2-3-3_0"
	echo "-p | --param: parameter file within the folder \"paramfiles/\". Default: default.txt"
	echo "-s | --seeds: seeds to run the simulator for each stage (e.g., 22). Default: 0."
  echo "-t | --timeout: timeout for each stage. Default: based on the number of nurses in the instance (60 + 6 * number of nurses * number of weeks)."
  echo "-o | --output: directory for the output. Default: outfiles/{param}/{instance}/{timestamp}."
	echo "-q | --queue: queue for sungrid. Default: empty."
}

# load arguments
instance_description="n030w4_1-2-3-3_0"
param="default.txt"
seeds=0
while [ ! -z $1 ]; do
  case $1 in
    -h|--help) printBashUsage
      exit 1;;
    -i | --instance) instance_description=$2; shift 2;;
		-p | --param) param=$2; shift 2;;
    -s | --seeds) seeds=$2; shift 2;;
    -t | --timeout) timeout=$2; shift 2;;
    -o | --output) outputDir=$2; shift 2;;
		-q | --queue) queue=$2; shift 2;;
    -*|--*) echo "Option unknown: $1"; shift 2;;
    *) echo "Cannot parse this argument: $1"
      printBashUsage
      exit 2;;
  esac
done

# parse the input competition instance name
parse=(`echo $instance_description | tr '_' ' ' `)
arrayArgs=(`echo ${parse[*]} | tr '-' ' ' `)
instance=${arrayArgs[0]}
numHist=${arrayArgs[1]}
nbArgs=${#arrayArgs[@]}
weekList=${parse[2]}
nbWeeks=$((nbArgs - 2))
if test ! -d "datasets/INRC2/$instance" ; then
	echo "Directory \"datasets/INRC2/$instance\" doest nos exist."
	exit 2;
fi

# get the parameters file from the second input
paramFile="paramfiles/$param"
if test ! -f $paramFile ; then
	echo "File \"$paramFile\" doest nos exist."
	exit 2;
fi
parseParam=(`echo $param | tr '.' ' ' `)
paramName=${parseParam[0]}

# display the inputs
echo "Instance: $instance_description"
echo "History number: ${parse[1]}"
echo "Week list: ${parse[2]}"
echo "Parsed instance: ${arrayArgs[*]}"
echo "Parameter file: ${paramFile}"

# create the output directory if it does not exist
if [ -z $outputDir ]; then
	timestamp=$(date +%s)
	outputDir="outfiles/${paramName}/$instance/$timestamp"
fi
mkdir -p $outputDir

# create the files that are going to be used in the execution
scenarioFile="datasets/INRC2/${instance}/Sc-${instance}.txt"
historyFile="datasets/INRC2/${instance}/H0-${instance}-${numHist}.txt"

for ((i=2; i<$nbArgs; i++)); do
	demandFiles[$i-2]="datasets/INRC2/${instance}/WD-${instance}-${arrayArgs[$i]}.txt"
done

for ((i=0; i<$nbWeeks; i++)); do
	solutionFiles[$i]="${outputDir}/sol-week${i}.txt"
done
validatorLog="${outputDir}/validator.txt"

# set the timeout depending on the operating system and number of nurses
if [ -z $timeout ]; then
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
fi
if  [ "$paramName" == "optimality" ] ; then
	timeout=86400
	echo "WARNING: Timeout changed to 86400 as parameters optimality used."
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
mkdir -p "bashfiles/${paramName}"
bashfile="bashfiles/${paramName}/$instance_description.sh"

# write the instructions in the bashfile
echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q $queue
#
# optimal script: run static scheduler" > ${bashfile}

echo "./bin/staticscheduler --dir datasets/INRC2/ --instance ${instance} --weeks ${weekList} --his ${numHist} --param ${paramFile} --sol ${outputDir} --timeout ${timeout} 2> ${outputDir}/log.txt"  >> ${bashfile}
echo "java -jar validator.jar --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${solutionFiles[*]} > ${validatorLog}" >> ${bashfile}

chmod 755 "${bashfile}"
