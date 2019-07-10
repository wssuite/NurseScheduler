#!/bin/bash

REPLICATION=30
FACTORS=( 1 2 4 6 8 10 )

INSTANCE="n100w4_0_1-1-0-8" # "n100w4_0_1-1-0-8" "n100w8_0_0-1-7-8-9-1-5-4"
TIME=2500
BASE_TIME=225 # 90% / 10 of TIME

BASE_OUTPUT="outfiles/${INSTANCE}/eval"
RESULTS="$BASE_OUTPUT/results.txt"

CURRENT_DIR="$(pwd)"
BASH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $BASH_DIR && cd ..

mkdir -p $BASE_OUTPUT
if [ ! -z $RESULTS ]; then
    echo "nEvaluationDemands time output validator_cost" >> $RESULTS
fi

for i1 in "${FACTORS[@]}"; do
  echo "i1: $i1"
  for i2 in "${FACTORS[@]}"; do
    GENERATION_TIME=$(( ${BASE_TIME} * $i2 ))
    echo "Generation time: $GENERATION_TIME, i2: $i2"
    for i3 in $(eval echo {1..$REPLICATION}); do
      OUTPUT="${BASE_OUTPUT}/${i1}_${GENERATION_TIME}_${i3}"
      echo "Create $OUTPUT"
      mkdir -p $OUTPUT

      echo "nEvaluationDemands=$i1" >> "$OUTPUT/stochasticOptions.txt"
      echo "solveToOptimality=1" >> "$OUTPUT/generationOptions.txt"
      echo "maxSolvingTimeSeconds=$GENERATION_TIME" >> "$OUTPUT/generationOptions.txt"

      # unset seeds from previous loop
      seeds=""
      source ./writeDynamicRun.sh -i $INSTANCE -t $TIME -o $OUTPUT

      echo "Run: ./${scriptfile}:"
      chmod +x *.jar
      cp validator.jar ./bin
      ./${scriptfile}
      rm -f ./${scriptfile}

      if [ ! -z $RESULTS ]; then
            COST=$(cat ${OUTPUT}/validator.txt | grep "Total cost:")
            echo "$i1 $TIME $OUTPUT $COST" >> $RESULTS
      fi
    done
  done
done


cd $CURRENT_DIR
