# !/bin/sh

# parse the input competition instance name
echo "Competition instance: $1"
parse=(`echo $1 | tr '_' ' ' `)
arrayArgs=(`echo ${parse[*]} | tr '-' ' ' `)
echo "${arrayArgs[*]}"
instance=${arrayArgs[0]}
numHist=${arrayArgs[1]}
nbArgs=${#arrayArgs[@]}
nbWeeks=$((nbArgs - 2))

# create the output directory if it does not exist
outputDir="outfiles/${1}/"
if test ! -d "outfiles" ; then
	echo "Create output directory"
	mkdir "outfiles"
fi
if test ! -d "${outputDir}" ; then
	echo "Create output directory"
	mkdir "${outputDir}"
fi
outputDir="outfiles/${1}/WithNoEvaluation/"
if test ! -d "${outputDir}" ; then
	echo "Create output directory"
	mkdir "${outputDir}"
fi

# create the files that are going to be used in the execution
scenarioFile="datasets/${instance}/Sc-${instance}.txt"
historyFile="datasets/${instance}/H0-${instance}-${numHist}.txt"


for ((i=0; i<=$nbWeeks; i++)); do
	solutionFiles[$i]="${outputDir}sol-week${i}.txt"
done
validatorLog="${outputDir}validatorOutput.txt"

for ((i=2; i<$nbArgs; i++)); do
	demandFiles[$i-2]="datasets/${instance}/WD-${instance}-${arrayArgs[$i]}.txt"
done

# print the files that are going to be used during the execution
echo "Instance: ${instance}"
echo "Number of weeks: ${nbWeeks}"
echo "Scenario file: ${scenarioFile}"
echo "History file: ${historyFile}"
echo "Week files: ${demandFiles[*]}"
echo "Output directory: ${outputDir}"
echo "Solution files: ${solutionFiles[*]}"
echo "Validator log file: ${validatorLog}"



sungridfile="bashfiles/sungridSimulator/${1}.sh"
chmod 755 "${sungridfile}"
echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator" > ${sungridfile}

echo "java -jar Simulator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --solver ./bin/roster --runDir ./ --outDir ${outputDir}" >> ${sungridfile}

#echo "java -jar validator.jar --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${solutionFiles[*]} > ${validatorLog}" >> ${sungridfile}
