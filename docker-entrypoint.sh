#!/bin/bash

# generate script
source generateScript.sh $1 $2 $3

# run script
chmod +x *.jar
cp validator.jar ./bin
./${scriptfile}
