#!/bin/bash

function printBashUsage {
  echo "This script will install all the depedencies."
  echo "Usage:"
  echo "-h | --help: display this message"
  echo "-b | --build-type: set the build type for cmake. Default: Relaease"
  echo "-d | --debug: set cmake build type to Debug."
  echo "-t | --tag: set docker image tag."
}

# load config arguments in one line
A=()
while [ ! -z "$1" ]; do
    for v in "$1"; do
        A+=("$v")
    done
    shift 1;
done

# parse arguments
BUILD="Release"
i=0
while [ ! -z ${A[${i}]} ]; do
  case ${A[${i}]} in
   -h | --help) printBashUsage
      exit 0;;
   -b | --build-type) BUILD="${A[((i+1))]}"; ((i+=2));;
   -d | --debug) BUILD="Debug"; ((i+=1));;
   -t | --tag) TAG=${A[((i+1))]}; ((i+=2));;
  esac
done

if [[ -z $TAG ]]; then
  if [[ "$BUILD" == "Release" ]]; then
    TAG="latest"
  else
    TAG="debug"
  fi
fi

# build docker image multiplatform
docker buildx build . -t legraina/nurse-scheduler:${TAG} --push \
       --platform linux/amd64,linux/arm64 --build-arg CMAKE_BUILD_TYPE=${BUILD}
