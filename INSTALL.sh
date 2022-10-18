#!/bin/bash

#CURRENT_DIR="$(pwd)"
#DIR=$1
#if [ -z $1 ]
#then
#  DIR=$CURRENT_DIR
#fi
#
## Install Boost
#wget -O boost_1_64_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.64.0/boost_1_64_0.tar.gz/download
#tar xzvf boost_1_64_0.tar.gz
#cd boost_1_64_0/
#./bootstrap.sh --exec-prefix=$DIR
#./b2 threading=multi
#./b2 install threading=multi
#cd ..
#rm -r boost_1_64_0
#rm boost_1_64_0.tar.gz

# Install BCP - look at the submodule
git submodule update --init --recursive
cd deps/coin
./install.sh