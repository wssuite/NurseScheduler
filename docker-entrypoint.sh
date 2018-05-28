#!/bin/bash

function printBashUsage {
  echo "This script will run the simulator and then the validator."
  echo "Usage:"
  echo "-h | --help: display this message"
  echo "-i | --instance: instance to simulate (must follow the pattern (data_weeks_history)). Default: n005w4_1-2-3-3_0"
  echo "-sc | --stochastic-config: config file for the stochastic parameters. Default: none"
  echo "-gc | --generation-config: config file for the generation parameters. Default: none"
  echo "-ec | --evaluation-config: config file for the evaluation parameters. Default: none"
  echo "-s | --seeds: seeds to run the simulator for each stage (e.g., 22-36-96-5). Default: random."
  echo "-t | --timeout: timeout for each stage. Default: based on the number of nurses in the instance."
  echo "-o | --output: directory for the output. Default. outfiles/{instance}/{seeds}_{timestamp}"
  echo "-g | --goal: goal to reach for the cost of the solution. Used for the unit tests. Default: none."
}

# generate script
source generateScript.sh "$@"

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
   -sc | --stochastic-config) cp "$2" "${outputDir}/stochasticOptions.txt"; shift 2;;
   -gc | --generation-config) cp "$2" "${outputDir}/generationOptions.txt"; shift 2;;
   -ec | --evaluation-config) cp "$2" "${outputDir}/evaluationOptions.txt"; shift 2;;
   -g | --goal) goal=$2; shift 2;;
   -*|--*) echo "Option unknown: $1"
      echo
      printBashUsage
      exit 2;;
   *) echo "Cannot parse this argument: $1"
      printBashUsage
      exit 2;;
  esac
done

# run script
echo "Run: ./${scriptfile}:"
cat "${scriptfile}"
chmod +x *.jar
cp validator.jar ./bin
./${scriptfile}

# if $3 defined, test the total cost
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
