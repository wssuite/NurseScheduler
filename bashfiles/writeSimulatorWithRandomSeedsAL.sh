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

# create the files that are going to be used in the execution
scenarioFile="datasets/${instance}/Sc-${instance}.txt"
historyFile="datasets/${instance}/H0-${instance}-${numHist}.txt"

for ((i=2; i<$nbArgs; i++)); do
	demandFiles[$i-2]="datasets/${instance}/WD-${instance}-${arrayArgs[$i]}.txt"
done

# create the root of output directory if it does not exist
outputDir="outfiles/Competition/${1}/"
if test ! -d "outfiles" ; then
	echo "Create output directory"
	mkdir "outfiles"
fi
if test ! -d "outfiles/Competition" ; then
        echo "Create output directory"
        mkdir "outfiles/Competition"
fi
if test ! -d "${outputDir}" ; then
	echo "Create output directory"
	mkdir "${outputDir}"
fi

# create the random seeds
for ((i=0; i<$nbWeeks; i++)); do
	# seeds[$i]=$RANDOM
	seeds[$i]=$(($RANDOM%100))
done

# create the specific output directory for these seeds
catseeds=${seeds[0]}
for ((i=1; i<$nbWeeks; i++)); do
	catseeds+="-${seeds[$i]}"
done

outputDir+="Flexible$catseeds"
if test ! -d "${outputDir}" ; then
	echo "Create output directory"
	mkdir "${outputDir}"
fi

sols="${outputDir}/sol-week0.txt"
for ((i=1; i<$nbWeeks; i++)); do
	sols+=" ${outputDir}/sol-week${i}.txt"
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
echo "solutions = ${sols}"
echo "log validator = ${outputDir}/validator.txt"

sungridfile="${2}/${1}_${catseeds}.sh"
echo "sungridfile=$sungridfile"
if test ! -d "${2}" ; then
	echo "Create run directory"
	mkdir "${2}"
fi

echo "#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -o /dev/null
#$ -q idra
#
# optimal script: launch the simulator
cd /home/legraina/git/RosterDesNurses/" > ${sungridfile}

echo "java -jar Simulator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --solver ./roster2 --runDir ./bin --outDir ${outputDir} --rand ${seeds[*]} --timeout ${timeout} --cus"  >> ${sungridfile}

echo "java -jar validator.jar  --sce ${scenarioFile} --his ${historyFile} --weeks ${demandFiles[*]} --sols ${sols} > ${outputDir}/validator.txt"  >> ${sungridfile}

chmod 755 "${sungridfile}"
