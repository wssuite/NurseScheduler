#!/bin/bash

function printBashUsage {
  echo "This script will run the simulator and then the validator."
  echo "Usage:"
  echo "-h | --help: display this message"
  echo "-d | --dynamic: Add this flag if you'd like to run the dynamic version."
  echo "-i | --instance: instance to simulate (must follow the pattern (data_weeks_history)). Default: n005w4_1-2-3-3_0"
  echo "-sc | --solver-config: config file for the solver parameters. Default: default or none if dynamic"
  echo "-gc | --generation-config: config file for the generation parameters. Default: none"
  echo "-ec | --evaluation-config: config file for the evaluation parameters. Default: none"
  echo "-s | --seeds: seeds to run the simulator for each stage (e.g., 22-36-96-5). Default: random."
  echo "-t | --timeout: timeout for the solver or the stage. Default: based on the number of nurses and weeks in the instance."
  echo "-o | --output: directory for the output. Default: outfiles/{instance}/{timestamp} or outfiles/{instance}/{seeds}_{timestamp} if dynamic"
  echo "-g | --goal: goal to reach for the cost of the solution. Used for the unit tests. Default: none."
}

# load config arguments
echo "$@"
while [ ! -z $1 ]; do
  case $1 in
    -h|--help) printBashUsage
      exit 0;;
   -i | --instance) instance_description=$2; shift 2;;
   -s | --seeds) seeds=$2; shift 2;;
   -t | --timeout) timeout=$2; shift 2;;
    # add config files
   -sc | --solver-config) param=$2; shift 2;;
   -gc | --generation-config) genParam=$2; shift 2;;
   -ec | --evaluation-config) evalParam=$2; shift 2;;
   -g | --goal) goal=$2; shift 2;;
   -d | --dynamic) dynamic="1"; shift 1;;
   -*|--*) echo "Option unknown: $1"
      echo
      printBashUsage
      exit 2;;
   *) echo "Cannot parse this argument: $1"
      printBashUsage
      exit 2;;
  esac
done

if [ -z $dynamic ]; then
	# parse inputs
	# parse the input competition instance name
	echo "Instance: $instance_description"
	parse=(`echo $instance_description | tr '_' ' ' `)
	instance=${parse[0]}
	hist=${parse[1]}
	weeks=${parse[2]}

	# create the root of output directory if it does not exist
	printf -v currenttime '%(%s)T' -1
	outputDir="outfiles/$instance_description/$currenttime"
	echo "Create output directory: ${outputDir}"
	mkdir -p "${outputDir}"

	# set default timeout
	if [ -z $timeout ]; then
		timeout=30
	fi

	# set default param file
	if [ -z param ]; then
		param="paramfiles/default.txt"
	fi

	# run the scheduler
	echo "Run: ./bin/staticscheduler --dir datasets/ --instance $instance --weeks $weeks --his $hist --param $param --sol $outputDir --timeout $timeout"
	./bin/staticscheduler --dir datasets/ --instance $instance --weeks $weeks --his $hist --param paramfiles/$param --sol $outputDir --timeout $timeout

	# run the validator
	./validator.sh $instance $weeks $hist $outputDir
else
	# generate script
	source ./scripts/writeDynamicRun.sh "$@"

    # copy config files
	if [ ! -z $param ]; then
		cp "$param" "${outputDir}/solverOptions.txt"
	fi
	if [ ! -z $genParam ]; then
		cp "$genParam" "${outputDir}/generationOptions.txt"
	fi
	if [ ! -z $evalParam ]; then
		cp "$evalParam" "${outputDir}/evaluationOptions.txt"
	fi

	# run script
	echo "Run: ./${scriptfile}:"
	cat "${scriptfile}"
	chmod +x *.jar
	cp validator.jar ./bin
	./${scriptfile}
fi

# if a goal is defined, test the total cost
if [ -z "$goal" ]
then
	exit 0
fi

# fetch the total cost and check if is the right one
cat ${outputDir}/validator.txt
result=$(cat ${outputDir}/validator.txt | grep "Total cost: $goal")
if [ -z "$result" ]
then
  echo "error"
  exit 1
fi
exit 0
