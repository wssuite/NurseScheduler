#!/bin/bash

function printBashUsage {
  echo "This script will run the simulator and then the validator."
  echo "Usage:"
  echo "-h | --help: display this message"
  echo "-d | --dynamic: Add this flag if you'd like to run the dynamic version."
  echo "-i | --instance: instance to simulate (must follow the pattern (data_weeks_history)). Default: n005w4_0_1-2-3-3"
  echo "-p, -sc | --param, --solver-config: config file for the solver parameters. Default: none"
  echo "-gc | --generation-config: config file for the generation parameters. Default: none"
  echo "-ec | --evaluation-config: config file for the evaluation parameters. Default: none"
  echo "-s | --seeds: seeds to run the simulator for each stage (e.g., 22-36-96-5). Default: random."
  echo "-t | --timeout: timeout for the solver or the stage. Default: based on the number of nurses and weeks in the instance."
  echo "-o | --output: directory for the output. Default: outfiles/{instance}/{timestamp} or outfiles/{instance}/{seeds}_{timestamp} if dynamic"
  echo "-g | --goal: goal to reach for the cost of the solution. Used for the unit tests. Default: none."
  echo "-v | --valgrind: use valgrind to run the code. Default: false."
  echo "-e | --evaluate: use the validator to evaluate thee solution. Default: true."
}

# load config arguments
echo "$@"
valgrindCMD="valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes"
instance_description="n005w4_0_1-2-3-3"
eval="1"
dynamic_args=""
other_args=""
while [ ! -z $1 ]; do
  case $1 in
    -h|--help) printBashUsage
      exit 0;;
   -i | --instance) instance_description=$2; shift 2;;
   -s | --seeds) seeds=$2; dynamic_args="${dynamic_args} -s $2"; shift 2;;
   -t | --timeout) timeout=$2; dynamic_args="${dynamic_args} -t $2"; shift 2;;
    # add config files
   -p | --param) param=$2; shift 2;;
   -sc | --solver-config) param=$2; shift 2;;
   -gc | --generation-config) genParam=$2; shift 2;;
   -ec | --evaluation-config) evalParam=$2; shift 2;;
   -g | --goal) goal=$2; shift 2;;
   -d | --dynamic) dynamic="1"; shift 1;;
   -v | -valgrind) valgrind="1"; shift 1;;
   -e | --evaluate) eval=$2; shift 2;;
   -*|--*) echo "Option unknown: $1. It will be passed to the scheduler."
      other_args="${other_args} $1 $2"; shift 2;;
   *) echo "Cannot parse this argument: $1"
      printBashUsage
      exit 2;;
  esac
done
dynamic_args="${dynamic_args} -i ${instance_description}"

if [ -z ${dynamic} ]; then
	# parse inputs
	# parse the input competition instance name
	echo "Instance: ${instance_description}"
	parse=(`echo ${instance_description} | tr '_' ' ' `)
	instance=${parse[0]}
	hist=${parse[1]}
	weeks=${parse[2]}

	# create the root of output directory if it does not exist
	currenttime=$( date +%s )
	outputDir="outfiles/${instance_description}/${currenttime}"
	echo "Create output directory: ${outputDir}"
	mkdir -p "${outputDir}"

	# base output
	sCMD="./bin/staticscheduler --dir datasets/ --instance ${instance} --weeks ${weeks} --his ${hist} --sol ${outputDir}"

	# set default timeout
	if [ -z ${timeout} ]; then
		timeout=30
	fi
	sCMD="${sCMD} --timeout ${timeout}"

	# set param file
	if [ ! -z param ]; then
		# param="paramfiles/default.txt"
		sCMD="${sCMD} --param paramfiles/${param}"
	fi

	# add others args
	sCMD="${sCMD} ${other_args}"

  if [ -z ${valgrind} ]; then
  	# run the scheduler
  	echo "Run: ${sCMD}"
  	${sCMD} | tee ${outputDir}/output.txt
    ret=${PIPESTATUS[0]}

  	# run the validator
  	echo ${ret}
  	if [ ${ret} -eq 0 -a ${eval} -eq 1 ]; then
  	    ./validator.sh ${instance} ${weeks} ${hist} ${outputDir} --verbose
  	fi
  else
    # run the scheduler with valgrind
  	echo "Run: ${valgrindCMD} ${sCMD}"
  	${valgrindCMD} ${sCMD}

    exit 0
  fi
else
	# generate script
	source ./scripts/writeDynamicRun.sh ${dynamic_args} ${other_args}

  # copy config files
	if [ ! -z ${param} ]; then
		cp "${param}" "${outputDir}/solverOptions.txt"
	fi
	if [ ! -z ${genParam} ]; then
		cp "${genParam}" "${outputDir}/generationOptions.txt"
	fi
	if [ ! -z ${evalParam} ]; then
		cp "${evalParam}" "${outputDir}/evaluationOptions.txt"
	fi

	# run script
	echo "Run: ./${scriptfile}:"
	cat "${scriptfile}"
	chmod +x *.jar
	cp validator.jar ./bin
	./${scriptfile}
	ret=${PIPESTATUS[0]}
    rm -f ${scriptfile}
fi

# display the solution
if [ ${ret} -eq 0 -a ${eval} -eq 1 ]; then
    cat ${outputDir}/validator.txt
fi

# if a goal is defined, test the total cost
if [ -z "$goal" ]
then
	exit 0
fi

# fetch the total cost and check if is the right one (in the bounds)
bounds=(`echo ${goal} | tr '-' ' ' `)
lb=${bounds[0]}
ub=${bounds[${#bounds[@]}-1]}
echo "lb=$lb; ub=$ub"

# check if should return an error (lb=-1)
if [ ${lb} == "E" ]; then
    if [ ${ret} -eq 0 ]; then
        exit 1
    fi
    exit 0
fi

# check bounds
if [ ${lb} -gt ${ub} ]; then
  echo "error: lb > ub"
  exit 1
fi

if [ ${eval} -eq 1 ]; then
    result=$(cat ${outputDir}/validator.txt | grep "Total cost:")
else
    result=$(cat ${outputDir}/output.txt | grep "Objective value")
fi
echo ${result}
if [ -z "$result" ]
then
  echo "error: total cost not found"
  exit 1
fi

rcost=$(echo ${result} | tr -dc '0-9')
echo "::set-output name=cost::${rcost}"
if [ ${rcost} -lt ${lb} ] || [ ${rcost} -gt ${ub} ]
then
  echo "error: bounds not respected"
  exit 1
fi

exit 0
