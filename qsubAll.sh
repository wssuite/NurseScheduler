#!/bin/bash

for entry in "sungrid"/*
do
  qsub ${entry}
done
