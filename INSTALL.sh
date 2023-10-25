#!/bin/bash

# Install BCP - look at the submodule
git submodule update --init --recursive
cd deps/coin
./install.sh