#!/bin/bash

valgrindOPT="--leak-check=full --show-leak-kinds=all --track-origins=yes"
dataDir="datasets/INRC2/"
instance_description="n005w4_0_1-2-3-3"
eval="1"
retries=0

function printBashUsage {
  echo "This script will run the simulator and then the validator."
  echo "Usage:"
  echo "-h  | --help: display this message"
  echo "-d  | --dynamic: Add this flag if you'd like to run the dynamic version."
  echo "-i  | --instance: instance to simulate (must follow the pattern (data_weeks_history)). Default: ${instance_description}"
  echo "-p  | --param: config file for the solver parameters. Default: none"
  echo "-sc | --solver-config: same than -p | --param"
  echo "-gc | --generation-config: config file for the generation parameters. Default: none"
  echo "-ec | --evaluation-config: config file for the evaluation parameters. Default: none"
  echo "-s  | --seeds: seeds to run the simulator for each stage (e.g., 22-36-96-5). Default: random."
  echo "-t  | --timeout: timeout for the solver or the stage. Default: based on the number of nurses and weeks in the instance."
  echo "-o  | --output: directory for the output. Default: outfiles/{instance}/{timestamp} or outfiles/{instance}/{seeds}_{timestamp} if dynamic"
  echo "-g  | --goal: goal to reach for the cost of the solution. Used for the unit tests. Default: none."
  echo "-v  | --valgrind: use valgrind to run the code. Default: false."
  echo "-vo | --valgrind-options: options for valgrind. Default: ${valgrindOPT}."
  echo "-e  | --evaluate: use the validator to evaluate thee solution. Default: ${eval}."
  echo "-r  | --root-dir-path: set the path where the script should be run. Default: do not move."
  echo "--retries: set the number of times to retry a failed run. Default: ${retries}."
  echo "--dir: set the dataset directory of the instance. Default: ${dataDir}."
}

# load config arguments in one line
ARGS=()
while [ ! -z "$1" ]; do
    for v in $1; do
        ARGS+=($v)
    done
    shift 1;
done
echo "${ARGS[@]}"
# parse arguments
dynamic_args=""
other_args=""
i=0
while [ ! -z ${ARGS[${i}]} ]; do
  case ${ARGS[${i}]} in
    -h|--help) printBashUsage
      exit 0;;
   -i | --instance) instance_description=${ARGS[((i+1))]}; ((i+=2));;
   -s | --seeds) seeds=${ARGS[((i+1))]}; dynamic_args="${dynamic_args} -s ${ARGS[((i+1))]}"; ((i+=2));;
   -t | --timeout) timeout=${ARGS[((i+1))]}; dynamic_args="${dynamic_args} -t ${ARGS[((i+1))]}"; ((i+=2));;
    # add config files
   -p | --param) param=${ARGS[((i+1))]}; ((i+=2));;
   -sc | --solver-config) param=${ARGS[((i+1))]}; ((i+=2));;
   -gc | --generation-config) genParam=${ARGS[((i+1))]}; ((i+=2));;
   -ec | --evaluation-config) evalParam=${ARGS[((i+1))]}; ((i+=2));;
   -g | --goal) goal=${ARGS[((i+1))]}; ((i+=2));;
   -d | --dynamic) dynamic="1"; ((i++));;
   -v | --valgrind) valgrind="1"; ((i++));;
   -vo | --valgrind-options) valgrindOPT=${ARGS[((i+1))]}; ((i+=2));;
   -e | --evaluate) eval=${ARGS[((i+1))]}; ((i+=2));;
   -r | --root-dir-path) rootDir=${ARGS[((i+1))]}; ((i+=2));;
   --pricer) pricer="1"; ((i+=1));;
   --retries) retries=${ARGS[((i+1))]}; ((i+=2));;
   --dir) dataDir=${ARGS[((i+1))]}; ((i+=2));;
   -*|--*) echo "Option unknown: ${ARGS[${i}]}. It will be passed to the scheduler."
      other_args="${other_args} ${ARGS[${i}]} ${ARGS[((i+1))]}"; ((i+=2));;
   *) echo "Cannot parse this argument: ${ARGS[${i}]}. It will be passed to the scheduler."
      other_args="${other_args} ${ARGS[${i}]}"; ((i+=1));;
  esac
done
dynamic_args="${dynamic_args} -i ${instance_description}"

# move to root dir if defined
if [ ! -z ${rootDir} ]; then
    cd ${rootDir}
fi

function run {
  # test pricer if defined
  if [ ! -z ${pricer} ]; then
    pCMD="./bin/pricer ${instance_description} ${other_args}"
    echo "Run: ${pCMD}"
    ${pCMD}
    return ${PIPESTATUS[0]}
  fi

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
    sCMD="./bin/staticscheduler --dir ${dataDir} --instance ${instance} --weeks ${weeks} --his ${hist} --sol ${outputDir}"

    # set default timeout
    if [ -z ${timeout} ]; then
      timeout=30
    fi
    sCMD="${sCMD} --timeout ${timeout}"

    # set param file
    if [ ! -z ${param} ]; then
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
      valgrindCMD="valgrind ${valgrindOPT}"
      echo "Run: ${valgrindCMD} ${sCMD}"
      ${valgrindCMD} ${sCMD}

      return 0
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
  if [ -z "$goal" ]; then
    return 0
  fi

  # fetch the total cost and check if is the right one (in the bounds)
  bounds=(`echo ${goal} | tr '-' ' ' `)
  lb=${bounds[0]}
  ub=${bounds[${#bounds[@]}-1]}
  echo "lb=$lb; ub=$ub"

  # check if should return an error (lb=-1)
  if [ ${lb} == "E" ]; then
      if [ ${ret} -eq 0 ]; then
          return 1
      fi
      return 0
  fi

  # check bounds
  if [ ${lb} -gt ${ub} ]; then
    echo "error: lb > ub"
    return 1
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
    return 1
  fi

  rcost=$(echo ${result} | tr -dc '0-9')
  echo "::set-output name=cost::${rcost}"
  if [ ${rcost} -lt ${lb} ] || [ ${rcost} -gt ${ub} ]
  then
    echo "error: bounds not respected"
    return 1
  fi

  return 0
}

ret=1
try=0
while [[ $try -le $retries && $ret -gt 0 ]]
do
  run
  ret=$?
  ((try++))
  echo "Try: $try and exit code: $ret"
done

exit $ret