#!/bin/bash
for bashfile in "bashfiles/sungridSimulator"/*
do
  qsub ${bashfile}
done
