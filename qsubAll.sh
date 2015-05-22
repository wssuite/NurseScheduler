#!/bin/bash

for entry in "sungrid"/*
do
  qsub ${entry} $1
done
