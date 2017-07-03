#!/bin/bash
for bashfile in ./bashfiles/${1}/*
do
	echo ${bashfile}
	${bashfile}
done
