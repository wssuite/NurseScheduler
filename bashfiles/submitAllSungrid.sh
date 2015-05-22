#!/bin/bash
cd ..
for bashfile in "bashfiles/sungridSimulator"/*
do
  qsub ${bahsfile}
done
