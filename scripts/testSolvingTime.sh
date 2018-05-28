#!/bin/bash

REPLICATION=30
FACTORS=10

INSTANCE="n100w4_0_1-1-0-8" # "n100w4_0_1-1-0-8" "n100w8_0_0-1-7-8-9-1-5-4"
nbNurses=100
BASE_TIME=$((($nbNurses-20)*3+10))

N_REL_DIVE=1
N_MIN_DIVE=2

BASE_OUTPUT="outfiles/${INSTANCE}/eval"
RESULTS="$BASE_OUTPUT/results.txt"

CURRENT_DIR="$(pwd)"
BASH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $BASH_DIR && cd ..

mkdir -p $BASE_OUTPUT
echo "nEvaluationDemands time output validator_cost" >> $RESULTS

for i1 in $(eval echo {1..$FACTORS}); do
  echo "i1: $i1"
  for i2 in $(eval echo {1..$FACTORS}); do
    TIME=$(($i2 * ${BASE_TIME}))
    echo "Time: $TIME, i2: $i2"
    for i3 in $(eval echo {1..$REPLICATION}); do
      OUTPUT="${BASE_OUTPUT}/${i1}_${TIME}_${i3}"
      echo "Create $OUTPUT"
      mkdir -p $OUTPUT

      echo "nEvaluationDemands=$i1" >> "$OUTPUT/stochasticOptions.txt"
      # echo "nbDiveIfRelGap=$(($i * $N_REL_DIVE))" >> "$OUTPUT/generationOptions.txt"
      # echo "nbDiveIfMinGap=$(($i * $N_ABS_DIVE))" >> "$OUTPUT/generationOptions.txt"

      source generateScript.sh -i $INSTANCE -t $TIME -o $OUTPUT

      echo "Run: ./${scriptfile}:"
      chmod +x *.jar
      cp validator.jar ./bin
      ./${scriptfile}

      COST=$(cat ${OUTPUT}/validator.txt | grep "Total cost:")
      echo "$i1 $TIME $OUTPUT $COST" >> $RESULTS
    done
  done
done


cd $CURRENT_DIR
