#!/bin/bash

CURRENT_DIR="$(pwd)"
DIR=$1
if [ -z $1 ]
then
  DIR=$CURRENT_DIR
fi

# Install BCP
cd $DIR
svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
mkdir -p Bcp-1.4/build
BCP_ARGS="--enable-shared=yes --enable-static=yes --with-cpx=no"
BUILD_TYPE="Release"
if [ "$2" == "Debug" ]; then
  BCP_ARGS="$BCP_ARGS --enable-debug"
fi
echo "Build type: $BUILD_TYPE; Run: ../configure $BCP_ARGS"
cd Bcp-1.4/build && ../configure $BCP_ARGS && make && make test && make install
cd $CURRENT_DIR
