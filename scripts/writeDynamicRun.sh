# !/bin/sh


function printBashUsage {
  echo "This script will run the simulator and then the validator."
  echo "Usage:"
  echo "-h | --help: display this message"
  echo "-i | --instance: instance to simulate (must follow the pattern (data_weeks_history)). Default: n030w4_1-2-3-3_0"
  echo "-s | --seeds: seeds to run the simulator for each stage (e.g., 22-36-96-5). Default: \$RANDOM."
  echo "-t | --timeout: timeout for each stage. Default: based on the number of nurses in the instance."
  echo "-o | --output: directory for the output. Default. outfiles/{instance}/{seeds}_{timestamp}"
}


# load arguments
instance_description="n030w4_1-2-3-3_0"
while [ ! -z $1 ]; do
  case $1 in
    -h|--help) printBashUsage
      exit 0;;
    -i | --instance) instance_description=$2; shift 2;;
    -s | --seeds) seeds=$2; shift 2;;
    -t | --timeout) timeout=$2; shift 2;;
    -o | --output) outputDir=$2; shift 2;;
    -*|--*) echo "Option unknown: $1"; shift 2;;
    *) echo "Cannot parse this argument: $1"
      printBashUsage
      exit 2;;
  esac
done

# parse the input competition instance name
echo "Competition instance: ${instance_description}"
parse=(`echo ${instance_description} | tr '_' ' ' `)
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

if [ -z ${seeds} ]; then
	# create the random seeds
	seeds=$(($RANDOM%100))
	for ((i=1; i<$nbWeeks; i++)); do
		seeds+="-$(($RANDOM%100))"
	done
fi
# get the seeds
allseeds=(`echo ${seeds} | tr '-' ' ' `)

# create the root of output directory if it does not exist
timestamp=$(date +%s)
if [ -z $outputDir ]; then
  outputDir="outfiles/${instance_description}/${seeds}_${timestamp}"
fi
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
if [ -z ${timeout} ]; then
    timeout=$((($nbNurses-20)*3+10))
fi
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
echo "seeds = ${allseeds[*]}"
echo "solutions = ${sols}"
echo "log validator = ${outputDir}/validator.txt"
echo "${allseeds[*]}" > "${outputDir}/seeds.txt"

# script file
scriptfile="${instance_description}_${seeds}_${timestamp}.sh"
echo "scriptfile=$scriptfile"

echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
" >> ${scriptfile}

echo "java -jar Simulator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --solver ./dynamicscheduler --runDir ./bin --outDir ${outputDir} --rand ${allseeds[*]} --timeout ${timeout} --cus"  >> ${scriptfile}

echo "java -jar validator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${sols} > ${outputDir}/validator.txt"  >> ${scriptfile}

chmod +x "${scriptfile}"
