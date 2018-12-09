#!/bin/bash

CURRENT_DIR="$(pwd)"
DIR=$1
if [ -z $1 ]
then
  DIR=$CURRENT_DIR
fi

# Install Boost
wget -O boost_1_64_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.64.0/boost_1_64_0.tar.gz/download
tar xzvf boost_1_64_0.tar.gz
cd boost_1_64_0/
./bootstrap.sh --exec-prefix=$DIR
./b2 threading=multi
./b2 install threading=multi
cd ..
rm -r boost_1_64_0
rm boost_1_64_0.tar.gz

# Install BCP
cd $DIR
svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
mkdir -p Bcp-1.4/build
BCP_ARGS="--enable-shared=yes --enable-static=yes --with-cpx=no"
BUILD_TYPE=$2
if [ -z $BUILD_TYPE ]; then
  BUILD_TYPE="Release"
elif [ "$BUILD_TYPE" == "Debug" ]; then
  BCP_ARGS="$BCP_ARGS --enable-debug"
fi
echo "Build type: $BUILD_TYPE; Run: ../configure $BCP_ARGS"
cd Bcp-1.4/build && ../configure $BCP_ARGS && make && make test && make install
cd $CURRENT_DIR
